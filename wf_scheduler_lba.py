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
from astropy.table import vstack,join
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
from astropy.coordinates import SkyCoord
### commensal extras
from vex import Vex
from astropy.coordinates import Angle
from collections import defaultdict

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
cat_path = sys.argv[2]
experiment = sys.argv[1]
print(experiment)
without_epoch = experiment[:-1]
vexfile = '{}.vex'.format(experiment)
if not os.path.exists(vexfile):
    os.system('wget ftp://ftp.atnf.csiro.au/pub/people/vlbi/{}/{}/{}.vex'.format(without_epoch, experiment, experiment))
vexfile=Vex('{}.vex'.format(experiment))
freq = vexfile.freq


##### Standard LBA inputs
FoV_query = 0.25
radius = 30*u.arcmin
do_plots = 'True'
npc = 100
do_targeted = 'True'
cat_type = 'csv'
filter_flux = 'True'
filter_value = 0.5 ### 5sigma noise level of LBA
phs_centre_fov = 58.24/60. ##smearing FoV, using the calculated bandwidth, phases within this distance will be combined
filter_overlap = 'True'
PB_plots = [25,32,100]
output_correlation_list = 'True'
phase_centre_format = 'difx'         #[csv|difx|sfxc]
prefix = '{}'.format(experiment)
filter_distance = 'False'
MSSC_additions='False'
MSSC_value = 1000
clip_phase_centres='True'
sortby = 'nearest'   # brightest | nearest
filter_exclusion_rad = 'True'
exclusion_radius = 10./60 #These phase centers will be removed to protect PI source
######

############ locate fringe finders, phase calibrators, target sources
fringe_finders = []
phase_refs = []
Sources = []


##fringe finders will be next to each other
for i in range(len(vexfile.sched)-3): ### minus three to remove final scan which is two phase ref scans
    if vexfile.sched[i]['source']==vexfile.sched[i+1]['source'] and vexfile.sched[i]['scan'][0]['scan_sec']<300:
        fringe_finders.append(vexfile.sched[i]['source'])

fringe_finders = list(set(fringe_finders))


neighbors = defaultdict(set)
time = defaultdict(float)

### phase reference and sources are one after another
for i in range(len(vexfile.sched) - 1):
    s1 = vexfile.sched[i]['source']
    s2 = vexfile.sched[i+1]['source']
    t1 = vexfile.sched[i]['scan'][0]['scan_sec']
    t2 = vexfile.sched[i+1]['scan'][0]['scan_sec']

    time[s1]=t1
    time[s2]=t2

    if s1 != s2 and s1 not in fringe_finders and s2 not in fringe_finders:
        neighbors[s1].add(s2)
        neighbors[s2].add(s1)

    pairs = {tuple(sorted([s, list(neigh)[0]])) for s, neigh in neighbors.items() if len(neigh) == 1}

    pairs_with_times = {}

    for a, b in pairs:
        pairs_with_times[(a, b)] = (time[a], time[b])

items = list(pairs_with_times.items())

##source will be longer observed than phase reference sources
for i in range(len(items)):
    sources = items[i][0]
    times = items[i][1]
    if times[0]>times[1]:
        phase_refs.append(sources[1])
        Sources.append(sources[0])
    else:
        phase_refs.append(sources[0])
        Sources.append(sources[1])

print('Fringe Finders: ', fringe_finders)
print('Phase Refs: ', phase_refs)
print('Sources: ', Sources)


### now locate potential phase centers from RACs mid, EMU or VLASS
source_names = []
ras = []
decs = []
survey = []
number_phase_centers = []
number_unfiltered_pc = []

for i in range(len(Sources)):
    source_name =Sources[i]
    source_names.append(source_name)
    ra_center = Angle(vexfile.source[source_name]['ra']).degree
    dec_center = Angle(vexfile.source[source_name]['dec']).degree
    ras.append(ra_center)
    decs.append(dec_center)
    pointing_centre = [ra_center, dec_center]
    freq  = vexfile.freq

    #### RACS Mid
    casdatap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
    job = casdatap.launch_job_async("SELECT * FROM AS110.racs_mid_components_v01 where 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{},{},{}))".format(ra_center,dec_center,FoV_query))
    r = job.get_results()
    keep_cols = ['ra','dec','total_flux']
    r = r[keep_cols]
    RA_column = keep_cols[0]
    Dec_column = keep_cols[1]
    flux_column = keep_cols[2]
    catalogue='RACS_potential_phase_centers_{}.csv'.format(source_name)
    print("Number of potential phase centers in RACs Mid", len(r))

    #### EMU
    emu_cat = Table.read(cat_path + '/EMU_components1.fits')
    center = SkyCoord(ra_center*u.deg, dec_center*u.deg)
    emu_coords = SkyCoord(ra=emu_cat['ra_deg_cont'], dec=emu_cat['dec_deg_cont'])
    sep = emu_coords.separation(center)
    emu_potential_phase = emu_cat[np.where((sep<=0.25*u.deg))]
    print('Number of potential phase centers in EMU',len(emu_potential_phase))

    #### VLASS
    vlass_cat = Table.read(cat_path + '/CIRADA_VLASS2QLv2_table1_components.csv',format='csv')
    vlass_coords = SkyCoord(ra=vlass_cat['RA']*u.deg, dec=vlass_cat['DEC']*u.deg)
    vlass_sep = vlass_coords.separation(center)
    vlass_potential_phase = vlass_cat[np.where((vlass_sep<=0.25*u.deg))]
    print("Number of potential phase centers in VLASS", len(vlass_potential_phase))


    ### pick survey with most sources to proceed with, if equal go with EMU
    if len(emu_potential_phase) >= len(r) and len(emu_potential_phase) >= len(vlass_potential_phase):
       catalogue = emu_potential_phase
       surv = 'EMU'
       flux_column = 'flux_int'
       RA_column = 'ra_deg_cont'
       Dec_column = 'dec_deg_cont'
       survey.append('EMU')
       print('Using EMU')

    elif len(r) > len(emu_potential_phase) and len(r) > len(vlass_potential_phase):
       catalogue = r
       surv = 'RACs'
       survey.append('RACs')
       print('Using RACs')

    elif len(vlass_potential_phase) > len(r) and len(vlass_potential_phase) > len(emu_potential_phase):
       surv = 'VLASS'
       catalogue = vlass_potential_phase
       flux_column = 'Total_flux'
       RA_column = 'RA'
       Dec_column = 'DEC'
       survey.append('VLASS')
       print('Using VLASS')

    if do_targeted == 'True':
        ### Read in tables
        df = catalogue
        master_table = catalogue
        filtered_coordinates=[]
        if filter_flux == 'True':
            logging.info('Flux filtering. All sources above %.2e kept' % (filter_value))
            df = df[df[flux_column]>filter_value]
            logging.info('Flux filtered. Nphs reduced from %d to %d' % (len(master_table[RA_column]),len(df[RA_column])))
        coords = SkyCoord(df[RA_column],df[Dec_column],unit=('deg','deg'))
        fluxs = df[flux_column]
        if filter_distance == 'True':
            logging.info('Filtering by distance from phase centre. All sources further than %.1f\' from phase centre will be removed' % radius)
            pointing_centres = SkyCoord(pointing_centre[0],pointing_centre[1],unit=('deg','deg'))
            truth_array = pointing_centres.separation(coords).to(u.arcmin).value < radius
            df = df[truth_array]
        if filter_exclusion_rad == 'True':
            logging.info('Filtering by distance from phase centre via exclusion radius. All sources within %.1f\' of phase centre will be removed as these are within the exclusion radius' % exclusion_radius)
            pointing_centres = SkyCoord(pointing_centre[0],pointing_centre[1],unit=('deg','deg'))
            truth_array_2 = pointing_centres.separation(coords).to(u.arcmin).value > exclusion_radius
            df = df[truth_array_2]
            logging.info('Removed sources within exclusion radius and outside FoV. Nphs reduced from %d to %d' % (len(master_table[RA_column]),len(df[RA_column])))
            if MSSC_additions=='True':
                logging.info('Adding in bright sources (above %.1f) in prior catalogue.' % MSSC_value)
                truth_array_3 = (pointing_centres.separation(coords).to(u.arcmin).value > radius)&(df[flux_column]>MSSC_value)
                df = df[truth_array_3]
            logging.info('Distance filtered. Nphs reduced from %d to %d' % (len(master_table[RA_column]),len(df[RA_column])))
            coords = SkyCoord(df[RA_column],df[Dec_column],unit=('deg','deg'))
            fluxs = df[flux_column]
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
            fluxs = df[flux_column]
        if filter_overlap == 'True':
            logging.info('Overlap filtering. Reducing number of phase centres if there are FoV overlaps')
            filtered_coordinates = filter_table(coords,phs_centre_fov) ## Filter the coordinates
        flux = df[flux_column]
        df = build_filtered_table(coords,flux,filter=filter_overlap,filter_indices=filtered_coordinates)
        master_table.rename_columns([RA_column, Dec_column], ['RA','DEC'])
        df['RA']  = np.round(df['RA'], 5)
        df['DEC'] = np.round(df['DEC'], 5)
        master_table['RA']  = np.round(master_table['RA'], 5)
        master_table['DEC'] = np.round(master_table['DEC'], 5)
        if filter_overlap == 'True':
            logging.info('Overlap filtered. Nphs reduced from %d to %d' % (len(master_table['RA']),len(df['RA'])))
        if clip_phase_centres == 'True':
            if ((len(df) < npc)&(len(df2) >= npc)):
                logging.info('Adding %d extra sources as overlap filter reduced phs centers'%(npc-len(df)))
                df = vstack([df,df2[npc:npc+(npc-len(df))]])
        number_phase_centers.append(len(df))
        number_unfiltered_pc.append(len(master_table))
        df.write('{}_{}_confirmed_phase_centers.csv'.format(source_name, surv),format='csv',overwrite=True)


    if do_plots == 'True':
        logging.info('Plotting phase centres')
        centre_coords = [np.average(df['RA']),np.average(df['DEC'])]
        pixels = 5000.
        large_range = np.max([np.max(master_table['RA'])-np.min(master_table['RA']),np.max(master_table['DEC'])-np.min(master_table['DEC'])])*0.5
        w = generate_central_wcs(centre_coords,[large_range/pixels,large_range/pixels],[0,0])
        fig = plt.figure(1,figsize=(9,9))
        ax = fig.add_subplot(111, projection=w)
        ax.scatter(df['RA'],df['DEC'],c='k',marker='+',transform=ax.get_transform('world'),s=20,label='Phase centres')
        ax.scatter(master_table['RA'],master_table['DEC'],transform=ax.get_transform('world'),s=2,label='Source positions')
        leg1 = ax.legend(loc='upper left', bbox_to_anchor=(1.01, 0.6))
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
        plt.close()

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
T['Survey'] = survey
T['number_of_og_PC'] = number_unfiltered_pc
T['number_phase_centers'] = number_phase_centers
T.write('{}_number_of_phase_center_info.csv'.format(experiment),format='csv', overwrite=True)
