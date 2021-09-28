#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
import os, ast, sys
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
import matplotlib.pyplot as plt
from astropy import wcs
import astropy.units as u
from astropy.io import fits
from matplotlib.colors import SymLogNorm
import pandas as pd
from simulator_functions import *

inputs = headless('simulator_inputs.txt') ## replace with sys.argv


# In[12]:


pointing_centre = ast.literal_eval(inputs['field_centre'])
size = ast.literal_eval(inputs['mosaic_area'])
pb_fwhm = np.load('%s/PB_fit.npy'%(inputs['output_path'])) 
c = SkyCoord(pointing_centre[0], pointing_centre[1],unit=('hour','deg'))


# In[15]:


coords = mosaic_pointings_square(centre_ra=c.ra.deg,centre_dec=c.dec.deg,ra_fov=size[0]/np.cos(c.dec.rad), dec_fov=size[1],theta=0, pointing_file_path='%s/mosaic.csv'%(inputs['output_path']),pb_fwhm=pb_fwhm,spacing=float(inputs['mosaic_filling_factor']))


# In[16]:


df= pd.read_csv('%s/mosaic.csv'%(inputs['output_path']), delim_whitespace=True)


# In[21]:


w = generate_central_wcs([c.ra.deg,c.dec.deg],[1/60,1/60],[1,1])


# In[24]:


fig = plt.figure(1,figsize=(9,9))
ax = fig.add_subplot(111,projection=w)
ax.scatter(df['RA'],df['DEC'],c='k',s=2,transform=ax.get_transform('world'))
#ax.scatter(t['RA'],t['DEC'],transform=ax.get_transform('world'),s=1)


for i,j in enumerate(df.RANGE):
	if i == 0:
		label = r'Phase centre ($10\%$ smearing)'
	else:
		label = ''
	r = SphericalCircle((df.RA[i]*u.deg,df.DEC[i]*u.deg), j/2. * u.degree,
					 edgecolor='k', facecolor='none',ls='--',
					 transform=ax.get_transform('world'),label=label)
	ax.add_patch(r)
ax.set_ylabel('Declination (J2000)')
ax.set_xlabel('Right Ascension (J2000)')
fig.savefig('%s/mosaic_pointings.pdf'%inputs['output_path'])


# In[ ]:




