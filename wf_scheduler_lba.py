#!/usr/bin/env python3
### Defaults
import os, sys
import copy
### Plotter
import matplotlib.pyplot as plt
### Numerics
import numpy as np
import argparse
from datetime import datetime
startTime = datetime.now()
### Table stuff
from astropy.io import ascii
from astropy.table import vstack
import pandas as pd
### Coordinate stuff
from astropy.coordinates import SkyCoord
import astropy.units as u
import ast
from astropy.visualization.wcsaxes import SphericalCircle
from itertools import cycle
### For the emailer and logger
from wf_vlbi_functions import *
import logging
import scipy.constants as c
from matplotlib.lines import Line2D
import warnings
### For querying catalogue
from astropy.table import Table
import astropy.units as u
warnings.filterwarnings("ignore", module = "matplotlib" )
warnings.filterwarnings("ignore", category=u.UnitsWarning)
from astropy.coordinates import SkyCoord
### commensal extras
from vex import Vex
from astropy.coordinates import Angle
from lba_functions import *
import pprint


def set_inputs(inputs):
    '''Recommend changing to a .yaml input format and rewriting this mess'''
    global do_targeted, catalogue, cat_type, RA_column, Dec_column,\
            filter_flat_flux, filter_by_pb, filter_by_pb_nsigma, filter_value,\
            flux_column, phs_centre_fov, filter_overlap, do_plots, PB_plots,\
            freq, output_correlation_list, phase_centre_format,\
            pointing_centre, prefix, filter_distance, radius, MSSC_value,\
            MSSC_additions, npc, clip_phase_centres, sortby, exclusion_radius

    do_targeted = ast.literal_eval(inputs.get('do_targeted', 'True'))
    catalogue = str(inputs.get('catalogue', ""))
    cat_type = str(inputs.get('table_format', 'csv'))
    RA_column = str(inputs.get('RA_column', "RA"))
    Dec_column = str(inputs.get('Dec_column', "Dec"))
    filter_flat_flux = ast.literal_eval(inputs.get('filter_flat_flux', 'False'))
    filter_by_pb = ast.literal_eval(inputs.get('filter_by_pb', 'True'))
    filter_by_pb_nsigma = float(inputs.get('filter_by_pb_nsigma', 5.0))
    filter_value = float(inputs.get('filter_value', 0.5))
    flux_column = str(inputs.get('flux_column', 'total_flux'))
    phs_centre_fov = convert_frac_to_float(inputs.get('phs_centre_fov', '58.24/60.'))
    filter_overlap= ast.literal_eval(inputs.get('filter_overlap', 'True'))
    do_plots= ast.literal_eval(inputs.get('do_plots', 'True'))
    PB_plots = ast.literal_eval(inputs.get('PBs', '[12,22,64]'))
    freq = inputs.get('observing_frequency', -1)
    output_correlation_list = ast.literal_eval(inputs.get('write_correlation_list', 'True'))
    phase_centre_format = str(inputs.get('phase_centre_format', 'difx').split(','))
    pointing_centre = ast.literal_eval(inputs.get('pointing_centre', '[None,None]'))
    prefix = str(inputs.get('catalogue_prefix', 'test'))
    filter_distance = ast.literal_eval(inputs.get('filter_distance', 'False'))
    radius = float(inputs.get('radius', 20)) # arcmin
    MSSC_value = float(inputs.get('MSSC_flux', 1000))
    MSSC_additions=ast.literal_eval(inputs.get('MSSC_additions', 'False'))
    npc = int(inputs.get('nphasecentres', 80))
    clip_phase_centres = ast.literal_eval(inputs.get('clip_phase_centres', 'False'))
    sortby = str(inputs.get('sortby', 'nearest'))
    exclusion_radius = float(inputs.get('exclusion_radius', 10./60.))
    #exclusion_radius = float(inputs.get('exclusion_radius', 0))

parser = argparse.ArgumentParser(
        description="Compare Vex file and catalogues and"
        " select phase centres for widefield correlation.")
parser.add_argument("filenames", type=str, nargs='+',
                    help="'Input' file, or .vex file")
parser.add_argument(
        '--lba', '-l',
        dest='lba', action='store_true',
        help="Use default LBA Catalogues, requires .vex file and you must set $WFCAT"
                " to be a directory with catalogues")

args = parser.parse_args()

### Setup logger
log_name = "%s.log" % os.path.basename(__file__).split('.py')[0]
setup_logging_to_file(log_name)
FORMAT = '%(asctime)s'
logging.basicConfig(format=FORMAT)
root = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
root.addHandler(handler)
logging.info('Beginning %s' % os.path.basename(__file__))

#filename = args.filename
# Get inputs
inputfilename = None
vexfile_name = None
vexfile = None
if len(args.filenames) > 2:
    raise ValueError('Too many files provided - give one .vex and/or one .input')
for filename in args.filenames:
    if re.search(r'.vex$', filename):
        vexfile_name = filename
        vexfile = Vex(vexfile_name)
        logging.info(f'Using Vex {filename}')
    else:
        inputfilename = filename
        logging.info(f'Using Input file {filename}')

# lba catalogue selection requires some info from vex file
if args.lba and (vexfile is None):
    parser.error('--lba requires a .vex file')

# parse inputs
sources_use = ()
calibrators = set()
targets = set()
if vexfile is not None:
    # use default inputs that usually work well with vexfile only
    inputs = {}
    # extract the source info from the vex file
    #fringe_finders, phase_refs, targets = locate_sources(vexfile)
    targets, calibrators = locate_sources2(vexfile)
    sources_use = calibrators.union(targets)
    logging.info('Calibrator(s): %s', ', '.join(calibrators))
    logging.info('Target(s): %s', ', '.join(targets))
    logging.info('All Source(s): %s', ', '.join(sources_use))
else:
    sources_use = ['pointing_centre']
# get inputs from the input file if there is one, overriding default values
if inputfilename is not None:
    inputs = headless(inputfilename)
set_inputs(inputs)

# A few inputs need to be set from the vex file, if not set in the input
if vexfile is not None:
    # always take prefix from the vex file if one is given
    prefix = re.sub('.vex', '', vexfile_name)
    # need the observing frequency
    if freq < 0.0:
        freq = vexfile.freq
logging.info(f'prefix: {prefix}')

# some inputs don't have reasonable defaults and must be set:
if freq < 0.0:
    raise ValueError("The frequency is not set in input file or vex")

# what catalogues are we using?
df = {}
master_table = {}
if not args.lba:
    # source and catalogue come from the input file
    df = ascii.read(catalogue,format=cat_type)
    #logging.info(df.info())
    master_table = ascii.read(catalogue,format=cat_type)
else:
    # read the LBA catalogues if requested
    # will determine df and master_table for each source later
    cat_path = os.environ["WFCAT"]
    logging.info(f'Searching {cat_path} for LBA catalogues')
    lba_catalogues = get_lba_catalogues(cat_path)

logging.info('Source(s) to use: %s', ', '.join(sources_use))

# some variables to record info on selections
source_names = []
ras = []
decs = []
survey = []
number_phase_centers = []
number_unfiltered_pc = []

# For each source select a best catalogue and cross-match
for source_name in sources_use:
    source_names.append(source_name)
    if vexfile is not None:
        # get source coords from the vex file if there is one, otherwise it
        # comes from input file
        ra_center = Angle(vexfile.source[source_name]['ra']).degree
        dec_center = Angle(vexfile.source[source_name]['dec']).degree
        pointing_centre = [ra_center, dec_center]
    ras.append(ra_center)
    decs.append(dec_center)

    surv = catalogue
    if args.lba:
        # search LBA catalogues for best option
        df, RA_column, Dec_column, flux_column, surv = select_lba_catalogue(
                source_name, lba_catalogues, radius, pointing_centre)
        logging.info(
                f'Chosen LBA catalogue: {surv}, {RA_column}, {Dec_column}, {flux_column}')
    #df = catalogue
    master_table = copy.deepcopy(df)
    survey.append(surv)

    if filter_flat_flux:
        logging.info('Flux filtering. All sources above %.2e kept' % (filter_value))
        df = df[df[flux_column]>filter_value]
        logging.info('Flux filtered. Nphs reduced from %d to %d' 
                     % (len(master_table[RA_column]),len(df[RA_column])))

    if filter_by_pb:
        #vexfile = str(inputs['vexfile'])
        # WARNING: need a better way to determine the flux units for generic catalogues
        flux_unit = u.Unit(str(inputs.get('flux_unit', 'mJy')))
        logging.info('PB sensitivity filtering using %s. Removing sources whose '
                     'primary-beam-attenuated flux is below %.1f sigma'
                     % (os.path.basename(vexfile_name), filter_by_pb_nsigma))
        pb_coords = SkyCoord(df[RA_column], df[Dec_column], unit=('deg', 'deg'))
        pointing_centres = SkyCoord(pointing_centre[0], pointing_centre[1], unit=('deg', 'deg'))
        offsets = pointing_centres.separation(pb_coords).to(u.deg).value
        # Effective image rms at each source's offset = central rms / PB power.
        # expected_rms_from_vex divides by the primary-beam power for offset>0,
        # so comparing the (un-attenuated) source flux to nsigma*eff_rms is
        # equivalent to requiring PB-attenuated_flux > nsigma * central_rms.
        #eff_rms = expected_rms_from_vex(vexfile_name, frequency=freq, offset=offsets)  # Jy/beam
        eff_rms = expected_rms_from_vex(
                vexfile_name, frequency=freq, offset=offsets,
                source=source_name, mk5clip=True)  # Jy/beam
        threshold = filter_by_pb_nsigma * eff_rms                                 # Jy/beam
        logging.info('Estimated central rms for %s: %.7f mJy/beam'%(source_name,np.min(eff_rms)*1e3))
        #logging.info(f'PB sensitivity filter: {dict(zip(offsets, threshold*1e3))} mJy/beam')
        flux_jy = (np.asarray(df[flux_column], dtype=float) * flux_unit).to(u.Jy).value
        df = df[flux_jy > threshold]
        logging.info('PB sensitivity filtered. Nphs reduced from %d to %d'
                     % (len(master_table[RA_column]), len(df[RA_column])))
        #printvals = [f'{x:.3f}' for x in sorted(df[flux_column])[0:10]]
        #logging.info(f'Faintest 10 fluxes: {printvals}')
        with np.printoptions(suppress=True, precision=3):
            logging.info(f'Faintest 10 fluxes (mJy): {np.array(sorted(df[flux_column])[0:10])}')

    if filter_distance:
        logging.info(
                f'Filtering by distance from phase centre. All sources further'
                f' than {radius}\' from phase centre will be removed')
        pointing_centres = SkyCoord(
                pointing_centre[0], pointing_centre[1], unit=('deg','deg'))
        coords = SkyCoord(df[RA_column], df[Dec_column], unit=('deg','deg'))   ## Generate skycoord instance of fits file
        truth_array = pointing_centres.separation(coords).to(u.arcmin).value < radius
        if MSSC_additions:
            logging.info('Adding in bright sources (above %.1f) in prior catalogue.' 
                         % MSSC_value)
            truth_array_2 = (
                    pointing_centres.separation(coords).to(u.arcmin).value >
                    radius) & (df[flux_column]>MSSC_value
                    )
            truth_array[truth_array_2==True]=True
        df = df[truth_array]
        logging.info('Distance filtered. Nphs reduced from %d to %d' 
                     % (len(master_table[RA_column]),len(df[RA_column])))

    if exclusion_radius > 0 and source_name in targets:
        logging.info('Filtering by distance from pointing centre via exclusion radius.'
                     ' All sources within %.1f\' of phase centre will be removed'
                     ' as these are within the exclusion radius' 
                     % exclusion_radius)
        pointing_centres = SkyCoord(
                pointing_centre[0], pointing_centre[1], unit=('deg','deg'))
        coords = SkyCoord(df[RA_column], df[Dec_column], unit=('deg','deg'))
        truth_array_3 = (pointing_centres.separation(coords).to(u.arcmin).value 
                         > exclusion_radius)
        df = df[truth_array_3]
        logging.info('Removed sources within exclusion radius and outside FoV.'
                ' Nphs reduced from %d to %d' 
                % (len(master_table[RA_column]),len(df[RA_column]))
                )

    if clip_phase_centres:
        if sortby == 'brightest':
            df.sort(keys=flux_column, reverse=True)
            df2 = Table([df[RA_column],df[Dec_column]], names=('RA','DEC'))
            df = df[0:npc]
        elif sortby == 'nearest':
            pointing_centres = SkyCoord(
                    pointing_centre[0],pointing_centre[1],unit=('deg','deg'))
            coords = SkyCoord(
                    df[RA_column], df[Dec_column], unit=('deg','deg'))
            df['separation'] = pointing_centres.separation(coords).to(u.arcmin).value
            df.sort(keys='separation', reverse=False)
            df2 = Table([df[RA_column],df[Dec_column]], names=('RA','DEC'))
            df = df[0:npc]

    filtered_coordinates=[]
    if filter_overlap:
        logging.info('Overlap filtering.'
                ' Reducing number of phase centres if there are FoV overlaps')
        coords = SkyCoord(df[RA_column], df[Dec_column], unit=('deg','deg'))
        filtered_coordinates = filter_table(coords, phs_centre_fov) ## Filter the coordinates

    flux = df[flux_column]
    coords = SkyCoord(df[RA_column], df[Dec_column], unit=('deg','deg'))
    df = build_filtered_table(
            coords, flux, filter=filter_overlap,
            filter_indices=filtered_coordinates)
    #master_table.rename_columns([RA_column, Dec_column], ['RA','DEC'])
    #df['RA']  = np.round(df['RA'], 5)
    #df['DEC'] = np.round(df['DEC'], 5)
    #master_table['RA']  = np.round(master_table['RA'], 5)
    #master_table['DEC'] = np.round(master_table['DEC'], 5)

    if filter_overlap:
        logging.info(
                'Overlap filtered. Nphs reduced from %d to %d'
                 %(len(master_table[RA_column]),len(df['RA'])))

    if clip_phase_centres:
        if ((len(df) < npc)&(len(df2) >= npc)):
            logging.info(
                    'Adding %d extra sources as overlap filter reduced phscenters'
                    %(npc-len(df)))
            df = vstack([df,df2[npc:npc+(npc-len(df))]])

    number_phase_centers.append(len(df))
    number_unfiltered_pc.append(len(master_table))

    # jump ship if no phase centres (plotter doesn't like it).
    if len(df) == 0:
        logging.warning(f'{source_name} returned no phase centres!')
        continue

    if do_plots:
        logging.info('Plotting phase centres')
        #centre_coords = [np.average(df['RA']),np.average(df['DEC'])]
        centre_coords = [pointing_centre[0], pointing_centre[1]]
        #print(centre_coords)
        pixels = 5000.
        large_range = np.max(
                [np.max(master_table[RA_column]) 
                 - np.min(master_table[RA_column]),
                 np.max(master_table[Dec_column]) 
                 - np.min(master_table[Dec_column])]) * 0.5
        w = generate_central_wcs(
                centre_coords, [large_range/pixels,large_range/pixels], [0,0])
        fig = plt.figure(figsize=(9,9))
        ax = fig.add_subplot(111, projection=w)
        # mark the centre
        ax.plot(
                centre_coords[0], centre_coords[1], marker='+', color='red',
                markersize=10, markeredgewidth=1, transform=ax.get_transform('world'))
        ax.scatter(
                df['RA'], df['DEC'], c='k', marker='+',
                transform=ax.get_transform('world'), s=20, 
                label='Phase centres')
        #print(df['RA'],df['DEC'])
        ax.scatter(
                master_table[RA_column], master_table[Dec_column],
                transform=ax.get_transform('world'), s=2, 
                label=f'Source positions \n ({radius} arcmin)')
        leg1 = ax.legend(loc='upper left', bbox_to_anchor=(1.01, 0.6))
        #ax.plot(df['RA'],df['DEC'],'-',transform=ax.get_transform('world'))
        #ax.set_xlim(pixels/-2.,pixels/2.)
        #ax.set_ylim(pixels/-2.,pixels/2.)
        for i in range(len(df['RA'])):
            r = SphericalCircle((
                df['RA'][i] * u.deg, df['DEC'][i] * u.deg), 
                phs_centre_fov * u.arcmin, 
                edgecolor='k', facecolor='none',
                transform=ax.get_transform('world'))
            ax.add_patch(r)
        if len(PB_plots) != 0:
            lst = cycle(['r','b','k','g'])
            ls2 = cycle(['--','-.',':'])
            custom_lines = []
            handles = []
            for i,j in enumerate(PB_plots):
                iter1 = next(lst)
                iter2 = next(ls2)
                PB_fov = ((c.c/freq)/j)*(180./np.pi)
                r = SphericalCircle(
                        (float(pointing_centre[0])*u.deg,
                        float(pointing_centre[1])*u.deg), 
                        PB_fov/2. * u.degree, edgecolor=iter1, linestyle=iter2,
                        facecolor='none', lw=2,
                        transform=ax.get_transform('world'))
                ax.add_patch(r)
                custom_lines.append(Line2D([0], [0], color=iter1,ls=iter2, lw=4))
                handles.append(r'$%s\,$m'%j)
        legend1 = ax.legend(
                custom_lines, handles, loc='upper left', 
                bbox_to_anchor=(1.01, 0.45),
                title=r'Primary beam')
        fig.add_artist(leg1)
        #ax.add_artist(legend1)
        ax.coords[0].set_axislabel('Right Ascension (J2000)')
        ax.coords[1].set_axislabel('Declination (J2000)')
        #fig.savefig('%s/%s_correlation_plot.pdf'%(os.getcwd(),prefix),bbox_inches='tight')
        ax.set_title(f'{source_name}, {surv}')
        fig.savefig(
                '{}_{}_correlation_plot.pdf'.format(prefix,source_name),
                bbox_inches='tight')
        #plt.show()
    
    if output_correlation_list:
        # first sort on distance - useful for correlator
        pointing_centres = SkyCoord(
                pointing_centre[0], pointing_centre[1], unit=('deg','deg'))
        coords = SkyCoord(df['RA'], df['DEC'], unit=('deg','deg'))
        df['separation'] = pointing_centres.separation(coords).to(u.arcmin).value
        df.sort(keys='separation', reverse=False)

        df.write(
                '{}_{}_confirmed_phase_centers.csv'.format(source_name, surv),
                format='csv', overwrite=True)
        if 'csv' in phase_centre_format:
            logging.info('Writing %d phase centres into CSV format'%len(df))
            outfile = '{}_{}_correlation_params.csv'.format(prefix,source_name)
            ascii.write(
                    df, outfile,
                    format='csv', fast_writer=False,overwrite=True)
            logging.info(f'Complete... {outfile} has been written to the cwd')
        if 'sfxc' in phase_centre_format:
            logging.info('Writing %d phase centres into VEX format'%len(df))
            write_correlation_params(
                    prefix=prefix+'_'+source_name, table=df, correlator='sfxc',
                    source=source_name)
            logging.info('Complete... %s_%s_correlation_params.vex has been written to the cwd' 
                         % (source_name, prefix))
        if 'difx' in phase_centre_format:
            logging.info('Writing %d phase centres into V2D format'%len(df))
            write_correlation_params(
                    prefix=prefix+'_'+source_name, table=df, correlator='difx',
                    source_name=source_name)
            logging.info('Complete... %s_%s_correlation_params.v2d has been written to the cwd' 
                         % (source_name, prefix))

logging.info('All sources done!')
T = Table()
T['Source_Name'] = source_names
T['ra'] = ras
T['dec'] = decs
T['Survey'] = survey
T['number_of_og_PC'] = number_unfiltered_pc
T['number_phase_centers'] = number_phase_centers
T.write(
        '{}_number_of_phase_center_info.csv'.format(prefix), format='csv',
        overwrite=True)
