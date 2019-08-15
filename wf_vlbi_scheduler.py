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
import pandas as pd
### Coordinate stuff
from astropy.coordinates import SkyCoord
import astropy.units as u
import ast
from astropy.visualization.wcsaxes import SphericalCircle

### For the emailer and logger
from wf_vlbi_functions import *
import logging
### Setup logger
log_name = "%s.log" % os.path.basename(__file__).split('.py')[0]
setup_logging_to_file(log_name)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
logging.info('Beginning %s' % os.path.basename(__file__))


inputs = headless('WFVLBI_inputs.txt')
do_targeted = str(inputs['targeted'])
catalogue = str(inputs['catalog'])
cat_type = str(inputs['type'])
RA_column = str(inputs['RA_column'])
Dec_column = str(inputs['Dec_column'])
filter_flux = str(inputs['filter_flux'])
filter_value = float(inputs['filter_value'])
flux_column = str(inputs['flux_column'])
phs_centre_fov = convert_frac_to_float(inputs['phs_centre_fov'])
filter_overlap=str(inputs['filter_overlap'])
do_plots= str(inputs['do_plots'])
PB_plots= ast.literal_eval(inputs['PBs'])

if do_targeted == 'True':
    ### Read in tables
    df=ascii.read(catalogue,format=cat_type)
    master_table = ascii.read(catalogue,format=cat_type)
    filtered_coordinates=[]
    if filter_flux == 'True':
        logging.info('Flux filtering. All sources above %.2f kept' % (filter_value))
        df = df[df[flux_column]>filter_value]
    if filter_overlap == 'True':
        logging.info('Overlap filtering. Reducing number of phase centres if there are FoV overlaps')
        coords = SkyCoord(df[RA_column],df[Dec_column],unit=('deg','deg'))   ## Generate skycoord instance of fits file
        filtered_coordinates = filter_table(coords,phs_centre_fov) ## Filter the coordinates

    df = build_filtered_table(coords,filter=filter_overlap,filter_indices=filtered_coordinates, RA_col=RA_column,Dec_col=Dec_column)
    if filter_overlap == 'True':
        logging.info('Overlap filtered. Nphs reduced from %d to %d' % (len(master_table[RA_column]),len(df['RA'])))


if do_plots == 'True':
    centre_coords = [np.average(df['RA']),np.average(df['DEC'])]
    pixels = 5000.
    large_range = np.max([np.max(master_table[RA_column])-np.min(master_table[RA_column]),np.max(master_table[Dec_column])-np.min(master_table[Dec_column])])*1.5
    w = generate_central_wcs(centre_coords,[large_range/pixels,large_range/pixels],[0,0])
    fig = plt.figure(1,figsize=(9,9))
    ax = fig.add_subplot(111, projection=w)
    ax.scatter(df['RA'],df['DEC'],transform=ax.get_transform('world'))
    print(len(master_table))
    ax.scatter(master_table[RA_column],master_table[Dec_column],transform=ax.get_transform('world'),s=1)
    ax.plot(df['RA'],df['DEC'],'.',transform=ax.get_transform('world'))
    ax.set_xlim(pixels/-2.,pixels/2.)
    ax.set_ylim(pixels/-2.,pixels/2.)
    for i in range(len(df['RA'])):
        r = SphericalCircle((df['RA'][i] * u.deg, df['DEC'][i] * u.deg), phs_centre_fov * u.arcmin,
                     edgecolor='k', facecolor='none',
                     transform=ax.get_transform('world'))
        ax.add_patch(r)
    fig.savefig('test.pdf')
    #plt.show()
