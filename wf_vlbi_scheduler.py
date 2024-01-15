### Defaults
import os, re, sys
### Plotter
import matplotlib.pyplot as plt
### Numerics
import numpy as np
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
warnings.filterwarnings("ignore", module = "matplotlib" )

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


inputs = headless(sys.argv[1])
do_targeted = str(inputs['targeted'])
catalogue = str(inputs['catalog'])
cat_type = str(inputs['table_format'])
RA_column = str(inputs['RA_column'])
Dec_column = str(inputs['Dec_column'])
filter_flux = str(inputs['filter_flux'])
filter_value = float(inputs['filter_value'])
flux_column = str(inputs['flux_column'])
phs_centre_fov = convert_frac_to_float(inputs['phs_centre_fov'])
filter_overlap=str(inputs['filter_overlap'])
do_plots= str(inputs['do_plots'])
PB_plots= ast.literal_eval(inputs['PBs'])
freq = float(inputs['observing_frequency'])
output_correlation_list = str(inputs['write_correlation_list'])
phase_centre_format = str(inputs['phase_centre_format']).split(',')
pointing_centre = ast.literal_eval(inputs['pointing_centre'])
prefix = str(inputs['catalogue_prefix'])
filter_distance = str(inputs['filter_distance'])
radius = float(inputs['radius'])
MSSC_value = float(inputs['MSSC_flux'])
MSSC_additions=str(inputs['MSSC_additions'])
npc = int(inputs['nphasecentres'])

if do_targeted == 'True':
    ### Read in tables
    df=ascii.read(catalogue,format=cat_type)
    #logging.info(df.info())
    master_table = ascii.read(catalogue,format=cat_type)
    filtered_coordinates=[]
    if filter_flux == 'True':
        logging.info('Flux filtering. All sources above %.2e kept' % (filter_value))
        df = df[df[flux_column]>filter_value]
        logging.info('Flux filtered. Nphs reduced from %d to %d' % (len(master_table[RA_column]),len(df['RA'])))
    coords = SkyCoord(df[RA_column],df[Dec_column],unit=('deg','deg'))   ## Generate skycoord instance of fits file
    if filter_distance == 'True':
        logging.info('Filtering by distance from phase centre. All sources further than %.1f\' from phase centre will be removed' % radius)
        pointing_centres = SkyCoord(pointing_centre[0],pointing_centre[1],unit=('deg','deg'))
        truth_array = pointing_centres.separation(coords).to(u.arcmin).value < radius
        if MSSC_additions=='True':
            logging.info('Adding in bright sources (above %.1f) in prior catalogue.' % MSSC_value)
            truth_array_2 = (pointing_centres.separation(coords).to(u.arcmin).value > radius)&(df[flux_column]>MSSC_value)
            truth_array[truth_array_2==True]=True
        df = df[truth_array]
        logging.info('Distance filtered. Nphs reduced from %d to %d' % (len(master_table[RA_column]),len(df['RA'])))
        coords = SkyCoord(df[RA_column],df[Dec_column],unit=('deg','deg'))
    if inputs['clip_phase_centres'] == 'True':
        if inputs['sortby'] == 'brightest':
            df.sort(keys=flux_column,reverse=True)
            df2 = Table([df[RA_column],df[Dec_column]], names=('RA','DEC'))
            df = df[0:npc]
        if inputs['sortby'] == 'nearest':
            pointing_centres = SkyCoord(pointing_centre[0],pointing_centre[1],unit=('deg','deg'))
            df['separation'] = pointing_centres.separation(coords).to(u.arcmin).value
            df.sort(keys='separation',reverse=False)
            df2 = Table([df[RA_column],df[Dec_column]], names=('RA','DEC'))
            df = df[0:npc]
        coords = SkyCoord(df[RA_column],df[Dec_column],unit=('deg','deg'))
    if filter_overlap == 'True':
        logging.info('Overlap filtering. Reducing number of phase centres if there are FoV overlaps')
        filtered_coordinates = filter_table(coords,phs_centre_fov) ## Filter the coordinates
    df = build_filtered_table(coords,filter=filter_overlap,filter_indices=filtered_coordinates)
    if filter_overlap == 'True':
        logging.info('Overlap filtered. Nphs reduced from %d to %d' % (len(master_table[RA_column]),len(df['RA'])))
    if inputs['clip_phase_centres'] == 'True':
        if ((len(df) < npc)&(len(df2) >= npc)):
            logging.info('Adding %d extra sources as overlap filter reduced phs centers'%(npc-len(df)))
            df = vstack([df,df2[npc:npc+(npc-len(df))]])



if do_plots == 'True':
    logging.info('Plotting phase centres')
    centre_coords = [np.average(df['RA']),np.average(df['DEC'])]
    pixels = 5000.
    large_range = np.max([np.max(master_table[RA_column])-np.min(master_table[RA_column]),np.max(master_table[Dec_column])-np.min(master_table[Dec_column])])*0.5
    w = generate_central_wcs(centre_coords,[large_range/pixels,large_range/pixels],[0,0])
    fig = plt.figure(1,figsize=(9,9))
    ax = fig.add_subplot(111, projection=w)
    ax.scatter(df['RA'],df['DEC'],c='k',marker='+',transform=ax.get_transform('world'),s=20,label='Phase centres')
    ax.scatter(master_table[RA_column],master_table[Dec_column],transform=ax.get_transform('world'),s=2,label='Source positions')
    leg1 = ax.legend(loc='upper left', bbox_to_anchor=(1.01, 0.6))
    #ax.plot(df['RA'],df['DEC'],'-',transform=ax.get_transform('world'))
    #ax.set_xlim(pixels/-2.,pixels/2.)
    #ax.set_ylim(pixels/-2.,pixels/2.)
    for i in range(len(df['RA'])):
        r = SphericalCircle((df['RA'][i] * u.deg, df['DEC'][i] * u.deg), phs_centre_fov * u.arcmin,
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
            print(PB_fov)
            r = SphericalCircle((float(pointing_centre[0])*u.deg, float(pointing_centre[1])*u.deg), PB_fov/2. * u.degree,
                     edgecolor=iter1,linestyle=iter2, facecolor='none',lw=2,
                     transform=ax.get_transform('world'))
            ax.add_patch(r)
            custom_lines.append(Line2D([0], [0], color=iter1,ls=iter2, lw=4))
            handles.append(r'%s\,m'%j)
    legend1 = ax.legend(custom_lines, handles, loc='upper left', bbox_to_anchor=(1.01, 0.45), title=r'\textbf{Primary beam}')
    ax.add_artist(leg1)
    ax.coords[0].set_axislabel('Right Ascension (J2000)')
    ax.coords[1].set_axislabel('Declination (J2000)')
    fig.savefig('%s/%s_correlation_plot.pdf'%(os.getcwd(),prefix),bbox_inches='tight')
    #plt.show()

if output_correlation_list == 'True':
    print(phase_centre_format)
    if 'csv' in phase_centre_format:
        logging.info('Writing %d phase centres into CSV format'%len(df))
        ascii.write(df, '%s_correlation_params.csv'%prefix, format='csv', fast_writer=False,overwrite=True)
        logging.info('Complete... %s_correlation_params.csv has been written to the cwd' % prefix)
    if 'sfxc' in phase_centre_format:
        logging.info('Writing %d phase centres into VEX format'%len(df))
        write_correlation_params(prefix=prefix,table=df,correlator='sfxc')
        logging.info('Complete... %s_correlation_params.vex has been written to the cwd' % prefix)
    if 'difx' in phase_centre_format:
        logging.info('Writing %d phase centres into V2D format'%len(df))
        write_correlation_params(prefix=prefix,table=df,correlator='difx')
        logging.info('Complete... %s_correlation_params.v2d has been written to the cwd' % prefix)
    
