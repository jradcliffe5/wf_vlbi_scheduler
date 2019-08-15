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

### For the emailer and logger
from wf_vlbi_functions import *
import logging
### Setup logger
log_name = "%s.log" % os.path.basename(__file__).split('.py')[0]
setup_logging_to_file(log_name)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
logging.info('Beginning %s' % os.path.basename(__file__))


inputs = headless('WFVLBI_inputs.txt')
do_targeted = bool(inputs['targeted'])
catalogue = str(inputs['catalog'])
cat_type = str(inputs['type'])
RA_column = str(inputs['RA_column'])
Dec_column = str(inputs['Dec_column'])
filter_flux = str(inputs['filter_flux'])
filter_value = float(inputs['filter_value'])
flux_column = str(inputs['flux_column'])
phs_centre_fov = convert_frac_to_float(inputs['phs_centre_fov'])
filter_overlap=str(inputs['filter_overlap'])


if do_targeted == True:
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

