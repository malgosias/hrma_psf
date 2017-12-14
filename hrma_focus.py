#!/usr/bin/env /proj/sot/ska3/flight/bin/python
#
# hrma_focus.py: extract data and plot hrma focus related plots
# author: t. isobe (tisobe@cfa.harvard.edu)
# last update: May 02, 2016
#
# last update: Nov 22, 2017, malgosia
#              arc4gl -> arc5gl
#              rewrote set_interval() to match date/time format for arc5gl
#

import numpy as np
import os
import sys
import re
from glob import glob
from shutil import rmtree
from datetime import datetime
import random
from astropy.io import fits
import time
from Ska.Shell import getenv, bash, Spawn
from Chandra.Time import DateTime
import Ska.arc5gl
from astropy.table import Table
from astropy.io import ascii
from scipy.optimize import curve_fit
import fileinput
import logging

LOG_FILENAME = f'hrma_focus.log'
logging.basicConfig(filename=LOG_FILENAME, filemode='a',
                    format='%(message)s', level=logging.DEBUG)

# HOUSE_KEEPING_FILE = '../Scripts/house_keeping/hrma_focus_house_keeping.txt'
HOUSE_KEEPING_FILE = 'hrma_focus_house_keeping.txt'
# Temporary file, current table is written to it before
# appending it to the house keeping file
SRC2_DATA_FILE = f'src2_data.tmp'

# HTML_DIR_PATH = '/data/mta4/www/DAILY/mta_src/'
HTML_DIR_PATH = '.'

# Filtering parameters
SNR_LIM = 6  # signal to noise ratio threshold for storage
SNR_LIM_PLT = 15 # signal to noise ratio threshold for psf fitting
RMAJ_HRC_LIM = 500  # source ellipse major axis threshold (pixels)
RMAJ_ACIS_LIM = 15  # source ellipse major axis threshold (pixels)
DEFOC_LIM = 0.01  # defocus
NUM_SRC_PER_OBSID = 1  # why it used to be 5? use obsids with at least # of sources

# Fit samples with > 5 sources
NUM_FIT_SRC = 5
COEFF1 = 2.07e-5
COEFF2 = 2.88e-5

# Sim-z offset limits for various instruments
SIMZ_MAX_ACIS = -150  # mm
SIMZ_MAX_ACIS_I = -210  # mm
SIMZ_MIN_HRC = 50  # mm (100 in the webpage)
SIMZ_MAX_HRC_I = 200  # mm

logging.info("\nI'll apply these filters ...")
logging.info(f'DEFOCUS threshold {DEFOC_LIM}')
logging.info(f'SNR threshold {SNR_LIM}')
logging.info(f'RMAJ threshold for ACIS {RMAJ_ACIS_LIM} pixels')
logging.info(f'RMAJ threshold for HRC {RMAJ_HRC_LIM} pixels')
logging.info(f'At least {NUM_SRC_PER_OBSID} src detection per obsid')

ascdsenv = getenv('source /home/ascds/.ascrc -r release; source /home/mta/bin/reset_param', shell='tcsh')

# !! Try to avoid the code below

# read directory list
# TODO: uncomment when running as mta

# path = '/data/mta/Script/Python_script2.7/house_keeping/dir_list'
# path = 'dir_list'

# f = open(path, 'r')
# items = [line.strip() for line in f.readlines()]
# f.close()

# mta_directories = {}

# for item in items:
#     tmp = re.split(':', item)
#     mta_directory_name  = tmp[1].strip()
#     mta_directory_location = tmp[0].strip()
#     mta_directories[mta_directory_name] = mta_directory_location

# append  pathes to private folders to a python directory
# sys.path.append(mta_directories['bin_dir'])
# sys.path.append(mta_directories['mta_dir'])

# couple of things needed

# TODO: uncomment bin_data and dare when running as mta
# bin_data = '/data/mta4/MTA/data/'

# dare = mcf.get_val('.dare', dir=bin_data, lst=1)
## hakama = mcf.get_val('.hakama', dir=bin_data, lst=1)


def download_data(start=None, stop=None):
    """
    Extract data to compute HRMA focus plots.

    :param start, stop: any DateTime compatible time formats
    :returns: fetches acis*evt2.fits.gz, hrc*evt2.fits.gz files
              from the archive and stores them in the current
              directory
    """

    # Remove directory 'param' and all files that it contains
    if os.path.isdir('param'):
        rmtree('param')
    os.makedirs('param')

    # check whether previous fits files are still around
    # if yes, remove them
    fits_files = glob('*fits*')
    for file_ in fits_files:
        os.remove(file_)

    # if time interval is not defined, set the relevant 1 month interval
    if start is None or stop is None:
        [start, stop] = set_interval()

    start = DateTime(start).date
    stop = DateTime(stop).date

    # extract acis and hrc evt2 files
    for inst in ('acis', 'hrc'):
        run_arc(inst, start, stop)


def set_interval():
    """
    Set time inteval:
     * if today's day is < 10, the inteval starts on the 1st day of
       the previous month and ends on the 1st day of the current month
     * if today's day is >= 10, the inteval starts on the 1st day of
       the current month and ends on the 1st of the next month

    :returns: start, stop times DateTime().fits format
    """
    # Today's date
    today = DateTime()
    year, month, day = today.year, today.mon, today.day

    if day < 10:
        start = f'{year}-{1:0>2}-01T00:00:00'.format(year, month - 1)
        stop = '{0:}-{1:0>2}-01T00:00:00'.format(year, month)
        if month == 1:
            start = '{}-12-01T00:00:00'.format(year - 1)
    else:
        start = '{0:}-{1:0>2}-01T00:00:00'.format(year, month)
        stop = '{0:}-{1:0>2}-01T00:00:00'.format(year, month + 1)
        if month == 12:
            stop = '{}-01-01T00:00:00'.format(year + 1)

    return (start, stop)


def run_arc(inst, start, stop):
    """
    Run arc5gl and extract evt2 data for instrument defined with `inst`.
    :param inst: instrument, `acis` or `hrc`
    :param start, stop: start and stop times, DateTime objects
    """

    arc5 = Ska.arc5gl.Arc5gl()

    # remember that mta's ~/.arc5gl_pwd does not correspond to use mta
    # run this as user malgosia

    arc5.sendline('operation=retrieve')
    arc5.sendline('dataset=flight')
    arc5.sendline(f'detector={inst}')
    arc5.sendline('level=2')
    arc5.sendline('filetype=evt2')
    arc5.sendline(f'tstart={start}')
    arc5.sendline(f'tstop={stop}')
    arc5.sendline('go')

    del arc5


def run_celldetect():
    """
    Run celldetect for each data set (*fits.gz, except of ERs) in the current
    directory. Creates *src2.fits files for each dataset.
    """

    # Take all *fits.gz files except ERs (er_files, obsids 5****, 6****)
    all_files = glob('*.fits.gz')
    er_files = glob('acisf[5|6]*.fits.gz') + glob('hrcf[5|6]*.fits.gz')
    fits_files = list(set(all_files) - set(er_files))

    bash("/usr/bin/env PERL5LIB=")

    for infile in fits_files:
        # Handle only none grating observations
        grating = read_header_value_from_file(infile, 'GRATING')

        if grating == 'NONE':
            outfile = infile.replace('evt', 'src')
            outfile = outfile.replace('.gz', '')

            # mode=h ??
            cmd_str = f'celldetect infile={infile} outfile={outfile} > /dev/null'

            # run celldetect
            try:
                bash(cmd_str, env=ascdsenv)
            except:
                logging.info(f'Celldetect failed on {infile} file')
                pass


def extract_data_from_src2():
    """
    Extracts information from the *src2.fits files. Reads headers and
    srclist data. Applies filtering criteria (defocus, snr, rmaj for
    HRC and ACIS, number of sources per obsid) and creates a text file
    with define with SRC2_DATA_FILE. This file contains sources that
    pass the filtering criteria.

    This file can be later read in as an astropy Table with columns:
    'OBSID', 'TSTART', 'TSTOP', 'SIMX', 'SIMZ', 'X', 'Y', 'SNR', 'RAVG',
    'RND', 'ROTANG', 'PSF', 'DIST', 'ANGDIST', 'INSTRUME', 'DIST_ARCSEC',
    'PSF_ARCSEC'.

    !! Creates source file for the current month - need to add this
    to the house keeping file with mission long data

    !! if cron job runs every Saturday then two consequtive runs
    may produce the same data - are they added twice to the house
    keeping log of mission long data?
    """

    print('')
    print(f' See {LOG_FILENAME} for processing stats')
    print('')

    dtype_ = [('X', '<f8'), ('Y', '<f8'), ('SNR', '<f4'), ('RMAJ', '<f4'),
              ('RMIN', '<f4'), ('ROTANG', '<f4'), ('PSFRATIO', '<f4'),
              ('NET_COUNTS', '<f4'), ('OBS_ID', '<U5'), ('TSTART', '<f8'),
              ('TSTOP', '<f8'), ('FP_TEMP', '<f8'), ('RA_NOM', '<f8'),
              ('DEC_NOM', '<f8'), ('ROLL_NOM', '<f8'), ('DEFOCUS', '<f8'),
              ('SIM_X', '<f8'), ('SIM_Y', '<f8'), ('SIM_Z', '<f8'),
              ('INSTRUME', '<U6'), ('RNDS', '<f4'), ('RAVG', '<f4'),
              ('PSF', '<f4'), ('DIST', '<f8'), ('PSF_ARCSEC', '<f4'),
              ('DIST_ARCSEC', '<f8'), ('ANGDIST', '<f8')]

    # Keys in srclist header that will be extracted
    # !! Do we need ROLL_NOM?
    header_keys = ('OBS_ID',  'TSTART', 'TSTOP', 'FP_TEMP',
                   'RA_NOM', 'DEC_NOM', 'ROLL_NOM', 'DEFOCUS',
                   'SIM_X', 'SIM_Y', 'SIM_Z')

    src2_files = glob('*src2.fits')

    total = 0
    total_flt = 0
    srcs_final = None

    for src2 in src2_files:

        logging.info(f'\nReading {src2} file')

        # Read the srclist data and header
        hdulist = fits.open(src2)
        srclist_data = hdulist[1].data
        srclist_header = hdulist[1].header
        hdulist.close()

        # Preliminary filtering: number of detected sources, defocus, simz
        if len(srclist_data) < NUM_SRC_PER_OBSID:
            logging.info(f'< {NUM_SRC_PER_OBSID} sources detected, {src2} dropped')
            continue

        defocus = read_header_value_from_header(srclist_header, 'DEFOCUS')
        if defocus is None:
            logging.info(f'DEFOCUS keyword not found, {src2} dropped')
            continue

        if np.abs(defocus) > DEFOC_LIM:
            logging.info(f'DEFOCUS {defocus}, above limit, {src2} dropped')
            continue

        simz = read_header_value_from_header(srclist_header, 'SIM_Z')
        if simz is None:
            logging.info(f'SIM_Z keyword not found, {src2} dropped')
            continue

        params = get_ref_params(simz)
        if np.any([par is None for par in params.values()]):
            logging.info(f"SIM_Z = {params['simz']}, {src2} dropped")
            continue

        total = total + len(srclist_data)

        # Store data as an astropy Table for easy filtering
        sources = {}
        sources['X'] = srclist_data['X']
        sources['Y'] = srclist_data['Y']
        sources['SNR'] = srclist_data['SNR']
        tmp = srclist_data['R'].T
        sources['RMAJ'] = tmp[0]
        sources['RMIN'] = tmp[1]
        sources['ROTANG'] = srclist_data['ROTANG']
        sources['PSFRATIO'] = srclist_data['PSFRATIO']
        sources['NET_COUNTS'] = srclist_data['NET_COUNTS']
        sources = Table(sources)

        # SNR and RMAJ filtering for this obsid
        flt1 = (sources['SNR'] > SNR_LIM) & (sources['RMAJ'] < params['rmaj_lim'])

        # Filter the shape, i.e. cases with PFSRATIO or RMIN equal zero
        flt2 = (sources['PSFRATIO'] > 0) & (sources['RMAJ'] > 0) & (sources['RMIN'] > 0)

        srcs_flt = sources[flt1 & flt2]

        n_srcs_flt = len(srcs_flt)
        if n_srcs_flt == 0:
            logging.info(f'0 sources passed SNR and shape filters, {src2} dropped')
            continue

        total_flt = total_flt + n_srcs_flt

        # Add columns with header info
        for key in header_keys:
            val = read_header_value_from_header(srclist_header, key)
            srcs_flt[key] = [val] * n_srcs_flt
            logging.info(f'{key} {val}')

        srcs_flt['INSTRUME'] = [params['instrume']] * n_srcs_flt

        # Derive extra columns (RNDS, ROTANG, RAVG, PSF, PSF_ARCSEC,
        # DIST, DIST_ARCSEC, ANGDIST) and add them to the Table. Use
        # sky coords (do not apply px2arcsec conversion) while computing
        # PSF, DIST, RAVG.

        srcs_flt['RNDS'] = srcs_flt['RMAJ'] / srcs_flt['RMIN']
        ravgs = 0.5 * (srcs_flt['RMAJ'] + srcs_flt['RMIN'])
        psfs = ravgs / srcs_flt['PSFRATIO']
        srcs_flt['RAVG'] = ravgs  # in pixels
        srcs_flt['PSF'] = psfs  # in pixels

        srcs_flt['ROTANG'] = srcs_flt['ROTANG'] * np.pi / 180  # rad

        dx = srcs_flt['X'] - params['xref']
        dy = srcs_flt['Y'] - params['yref']
        srcs_flt['DIST'] = np.sqrt(dx**2 + dy**2)  # in pixels

        # Add PSF and DIST in arcsec
        srcs_flt['PSF_ARCSEC'] = srcs_flt['PSF'] * params['px2arcsec']
        srcs_flt['DIST_ARCSEC'] = srcs_flt['DIST'] * params['px2arcsec']

        angdists = np.arctan(dy / dx)
        for angdist in angdists:
            while angdist < 0:
                angdist = angdist + np.pi
        srcs_flt['ANGDIST'] = angdists

        logging.info(f'{len(sources)} sources read')
        logging.info(f'{n_srcs_flt} passed filtering criteria')

        # Stack filtered sources as arrays
        srcs_flt_array = np.array(srcs_flt.as_array(), dtype=dtype_)
        if srcs_final is None:
            srcs_final = srcs_flt_array
        else:
            srcs_final = np.hstack((srcs_final, srcs_flt_array))

    logging.info(f'\n{total} sources read (total)')
    logging.info(f'{total_flt} sources passed filters')

    # Write the new source data to the temporary SRC2_DATA_FILE,
    # then append this file to the house keeping file, skip the first line
    # (with column names) if the house keeping file already exists
    if total_flt > 0:
        ascii.write(srcs_final, SRC2_DATA_FILE, overwrite=True)
        line_number = 1
        if os.path.isfile(HOUSE_KEEPING_FILE):
            line_number = 2
        cmd_str = f'tail -n +{line_number} {SRC2_DATA_FILE} >> {HOUSE_KEEPING_FILE}'
        bash(cmd_str)


def get_ref_params(simz):
    """
    Get reference parameters for each obsid.

    :param simz: SIM z location read from the fits file header
    :returns: dictionary with keys 'xref', 'yref', 'px2arcsec',
              'rmaj_lim', 'instrume' representing aimpoint
              coordinates, pixel scale to convert from pixels
              to arcsec, limit on the major axis of the ellipse
              and the instrument name).
    """
    if simz < SIMZ_MAX_ACIS:  # acis
        xref = 4096.5
        yref = 4096.5
        px2arcsec = 0.492  # arcsec/pix
        rmaj_lim = RMAJ_ACIS_LIM

        if simz < SIMZ_MAX_ACIS_I:  # acis-I
           # zoff = -233.6 - simz
           instrume = 'ACIS-I'
        else:  # acis-S
           # zoff = -190.1 - simz
           instrume = 'ACIS-S'

    elif simz > SIMZ_MIN_HRC:  # hrc
        px2arcsec = 0.13175  # arcsec/pix
        rmaj_lim = RMAJ_HRC_LIM

        if simz < SIMZ_MAX_HRC_I:  # hrc-I
            # zoff = 126.99 - simz
            xref = 16384.5
            yref = 16384.5  # 0.006429 mm/pixel
            instrume = 'HRC-I'
        else:  # hrc-S
            # zoff = 250.1 - simz
            xref = 32768.5
            yref = 32768.5  # 0.006429 mm/pixel
            instrume = 'HRC-S'

    else:
        return [None] * 6

    logging.info(f'Aimpoint ({xref}, {yref})')

    params = {'xref': xref, 'yref': yref, 'px2arcsec': px2arcsec,
              'rmaj_lim': rmaj_lim, 'instrume': instrume}

    return params


def read_header_value_from_file(fits_file, name):
    """
    Read fits file header value for a given parameter `name`
    :param fits_file: name of the fits file
    :param name: name of the header keyword
    :returns: val, keyword value or None if keyword not found
    """
    hdu = fits.open(fits_file)
    header = hdu[1].header
    hdu.close()

    try:
        val = header[name.lower()]
    except:
        val = None

    return val


def read_header_value_from_header(header, name):
    """
    Find keyword value in a fits file header
    :param header: fits file header
    :param name: name of the header keyword
    :returns: val, keyword value, or None if keyword not found
    """
    try:
        val = header[name.lower()]
    except:
        val = None

    return val


def fit_psf(infile=None):
    """
    Fit the PSF
    """

    if infile is None:
        infile = HOUSE_KEEPING_FILE

    sources = ascii.read(infile)

    # Make plots using sources with snr > SNR_LIM_PLT
    flt = sources['SNR'] > SNR_LIM_PLT
    sources = sources[flt]

    psfs = []

    for label in ('all', 'ACIS-I', 'ACIS-S', 'HRC-I', 'HRC-S'):

        flt = sources['INSTRUME'] == label
        srcs_flt = sources[flt]

        if len(srcs_flt) > NUM_FIT_SRC:
            psf = PSF(label)
            psf.psf_time_anal(srcs_flt, label)
            psfs.append(psf)

    return psfs



def psf_function_a0(x, a0):
    """
    Fit only for a0 (on-axis PSF)
    """
    return a0 + COEFF1 * x + COEFF2 * x**2



def psf_function(x, a0, a1, a2):
    """
    Fit for all coefficients
    """
    return a0 + a1 * x + a2 * x**2



class PSF():

    x = None
    y = None
    popt = None
    pcov = None
    on_axis = None
    on_axis_error = None

    def __init__(self, instrume, year=None):
        self.instrume = instrume
        self.year = year


    def psf_time_anal(sources, title):
        """
        Time analysis - refit only the current year, somehow remember or
        read from file coefficients for the previous years? They will change
        if we change the source detection algorithm?
        """

        if self.year is None:
            current_year = DateTime().year
        # years = np.arange(1999, current_year + 1)
        years = [current_year]

        for year in years:
            self.year = year
            flt1 = DateTime(sources['TSTART']).date >= np.str(year)
            flt2 = DateTime(sources['TSTART']).date < np.str(year + 1)
            flt_sources = sources[flt1 & flt2]
            if len(flt_sources) > NUM_FIT_SRC:
               flt_sources.sort('DIST_ARCSEC')
               x = flt_sources['DIST_ARCSEC']
               y = flt_sources['PSF_ARCSEC']
               self.popt, self.pcov = curve_fit(psf_function, x, y)
               popt, pcov = curve_fit(psf_function_a0, x, y)
               self.on_axis = popt[0]
               self.on_axis_error = pcov[0][0]
               self.x = x
               self.y = y


def generate_plots(infile=None):
    """
    Generate all the plots, move them to the HTML_DIR_PATH
    and clean the fits files (clean the fits files earlier!)
    """

    if infile is None:
        infile = HOUSE_KEEPING_FILE

    sources = ascii.read(infile)

    # Make plots using sources with snr > SNR_LIM_PLT
    flt = sources['SNR'] > SNR_LIM_PLT
    sources = sources[flt]

    for label in ('all', 'ACIS-I', 'ACIS-S', 'HRC-I', 'HRC-S'):

        flt = sources['INSTRUME'] == label
        srcs_flt = sources[flt]

        psfs = []

        if len(srcs_flt) > 0:
            psf = PSF(label)
            psf.psf_time_anal(srcs_flt, label)
            psfs.append(psf)
            # src_mon_plots(srcs_flt, title=label)

            # plots by dist
            # rbins = [(0, 30), (30, 60), (60, 120), (120, 180), (180, 240),
            #          (240, 300), (300, 420), (420, 600)]

            # for rbin in rbins:
            #     flt1 = srcs_flt['DIST'] >= rbin[0]
            #     flt2 = srcs_flt['DIST'] < rbin[1]
            #     r_sources = srcs_flt[flt1 & flt2]
            #     if len(r_sources) > 0:
            #         title = '{}_r{}-{}'.format(label, rbin[0], rbin[1])
            #         src_mon_plots(r_sources, title=title)


    # Move the plots and remove the fits files
    cmd_str = 'mv -f *.html *.gif ' + HTML_DIR_PATH + '/. 2>/dev/null'
    os.system(cmd)

    os.system('rm -f *.fits')

    # Make webpages based on psfs



def src_mon_plots(sources, title):
    return


def update_index_page():
    """
    Update the "Last Updated" date in the index.html file
    """

    today = DateTime().caldate
    year, month, day = today[:4], today[4:7], today[7:10]

    html_file = HTML_DIR_PATH + 'index.html'

    # e.g. Last Updated: Nov 04, 2019
    newline = 'Last Updated: ' + month + ' ' + day + ', ' + year

    for line in fileinput.input(html_file, inplace=True):
        if 'Last Updated' in line:
            print(newline)
        else:
            print(line)
    fileinput.close()



if __name__ == "__main__":

    if len(sys.argv) > 2:
        start = sys.argv[1].strip()
        stop  = sys.argv[2].strip()
    else:
        start = ''
        stop  = ''

    download_data(start, stop)

    run_celldetect()

    extract_data_from_src2()

    # generate_plots()

    # update_index_page()
