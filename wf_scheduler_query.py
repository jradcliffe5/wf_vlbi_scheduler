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
### For querying catalogue
from astropy.table import Table
from astroquery.utils.tap.core import TapPlus
import astropy.units as u
warnings.filterwarnings("ignore", module = "matplotlib" )
###
from vex import Vex
from astropy.coordinates import Angle

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

#### Download VEX

experiment = sys.argv[1]
print(experiment)
without_epoch = experiment[:-1]
print(without_epoch)
os.system('wget ftp://ftp.atnf.csiro.au/pub/people/vlbi/{}/{}/{}.vex'.format(without_epoch, experiment, experiment))


vexfile=Vex('{}.vex'.format(experiment))

##### Make these standard for all LBA observations
FoV_query = 0.25 
radius = 30*u.arcmin
do_plots = 'True'
npc = 100
do_targeted = 'True'
cat_type = 'csv'
filter_flux = 'True'
filter_value = 0.5 ### 5sigma noise level of LBA
phs_centre_fov = 20./60. ##smearing FoV, just standard from Jacks example
filter_overlap = 'True'
PB_plots = [25,32,100] ### again just using Jacks example
output_correlation_list = 'True'
phase_centre_format = 'difx'         #[csv|difx|sfxc]
prefix = '{}'.format(experiment)
filter_distance = 'False'
MSSC_additions='False'
MSSC_value = 1000
clip_phase_centres='True'
sortby = 'nearest'   # brightest | nearest
######

Sources = vexfile.source
number_of_sources = len(list(vexfile.source.keys()))
print(number_of_sources)

source_names = []
ras = []
decs = []
number_phase_centers = []

for i in range(number_of_sources):
    source_name = list(vexfile.source.keys())[i]
    source_names.append(source_name)
    #first_source_name = list(vexfile.source.keys())[0] ### just for a text, will maybe loop through this?
    #source = vexfile.source[first_source_name]
    ra_center = Angle(vexfile.source[source_name]['ra']).degree
    dec_center = Angle(vexfile.source[source_name]['dec']).degree ### need to convert 
    ras.append(ra_center)
    decs.append(dec_center)
    pointing_centre = [ra_center, dec_center]
    freq  = vexfile.freq
    print(pointing_centre)

    casdatap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
    job = casdatap.launch_job_async("SELECT * FROM AS110.racs_mid_components_v01 where 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{},{},{}))".format(ra_center,dec_center,FoV_query))
    r = job.get_results()
    keep_cols = ['ra','dec','total_flux']
    r = r[keep_cols]
    print(r)
    r.write('RACS_potential_phase_centers_{}.csv'.format(source_name),overwrite=True)
    RA_column = keep_cols[0]
    Dec_column = keep_cols[1]
    flux_column = keep_cols[2]
    catalogue='RACS_potential_phase_centers_{}.csv'.format(source_name)

    if do_targeted == 'True':
        ### Read in tables
        df=ascii.read(catalogue,format=cat_type)
        #logging.info(df.info())
        master_table = ascii.read(catalogue,format=cat_type)
        filtered_coordinates=[]
        if filter_flux == 'True':
            logging.info('Flux filtering. All sources above %.2e kept' % (filter_value))
            df = df[df[flux_column]>filter_value]
            logging.info('Flux filtered. Nphs reduced from %d to %d' % (len(master_table[RA_column]),len(df[RA_column])))
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
            logging.info('Distance filtered. Nphs reduced from %d to %d' % (len(master_table[RA_column]),len(df[RA_column])))
            coords = SkyCoord(df[RA_column],df[Dec_column],unit=('deg','deg'))
        if clip_phase_centres == 'True':
            if sortby == 'brightest':
                df.sort(keys=flux_column,reverse=True)
                df2 = Table([df[RA_column],df[Dec_column]], names=('RA','DEC'))
                df = df[0:npc]
            if sortby == 'nearest':
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
        if clip_phase_centres == 'True':
            if ((len(df) < npc)&(len(df2) >= npc)):
                logging.info('Adding %d extra sources as overlap filter reduced phs centers'%(npc-len(df)))
                df = vstack([df,df2[npc:npc+(npc-len(df))]])
        number_phase_centers.append(len(df))



    if do_plots == 'True':
        logging.info('Plotting phase centres')
        centre_coords = [np.average(df['RA']),np.average(df['DEC'])]
        print(centre_coords)
        pixels = 5000.
        large_range = np.max([np.max(master_table[RA_column])-np.min(master_table[RA_column]),np.max(master_table[Dec_column])-np.min(master_table[Dec_column])])*0.5
        w = generate_central_wcs(centre_coords,[large_range/pixels,large_range/pixels],[0,0])
        fig = plt.figure(1,figsize=(9,9))
        ax = fig.add_subplot(111, projection=w)
        ax.scatter(df['RA'],df['DEC'],c='k',marker='+',transform=ax.get_transform('world'),s=20,label='Phase centres')
        print(df['RA'],df['DEC'])
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
        fig.savefig('{}_{}_correlation_plot.pdf'.format(prefix,source_name),bbox_inches='tight')
        #plt.show()

    if output_correlation_list == 'True':
        if 'csv' in phase_centre_format:
            logging.info('Writing %d phase centres into CSV format'%len(df))
            ascii.write(df, '{}_{}_correlation_params.csv'.format(prefix,source_name), format='csv', fast_writer=False,overwrite=True)
            logging.info('Complete... %s_correlation_params.csv has been written to the cwd' % prefix)
        if 'sfxc' in phase_centre_format:
            logging.info('Writing %d phase centres into VEX format'%len(df))
            write_correlation_params(prefix=prefix,table=df,correlator='sfxc')
            logging.info('Complete... %s_correlation_params.vex has been written to the cwd' % prefix)
        if 'difx' in phase_centre_format:
            logging.info('Writing %d phase centres into V2D format'%len(df))
            write_correlation_params(prefix=prefix+'_'+source_name,table=df,correlator='difx')
            logging.info('Complete... %s_correlation_params.v2d has been written to the cwd' % prefix)
    
T = Table()
T['Source_Name'] = source_names
T['ra'] = ras
T['dec'] = decs
T['number_phase_centers'] = number_phase_centers
T.write('{}_number_of_phase_center_info.csv'.format(experiment),format='csv', overwrite=True)
