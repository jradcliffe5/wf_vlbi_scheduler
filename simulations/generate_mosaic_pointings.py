import matplotlib
matplotlib.use('Agg')

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


def generate_rectangle(x_fov,y_fov,RA_centre,Dec_centre,theta):
	fov = np.array([x_fov/2.,y_fov/2.])
	xy = np.array([fov,fov*np.array([-1,1]),fov*-1,fov*np.array([1,-1])]).T
	return translate(rotate(xy,theta*(np.pi/180.)),[RA_centre,Dec_centre]).T

def mosaic_pointings_square(centre_ra, centre_dec, ra_fov, dec_fov, theta, pointing_file_path, pb_fwhm, spacing=0.6):

	# Initialise some things...

	ra_mosaic_buffer = []   # array of RA coordinates for plotting   
	dec_mosaic_buffer = []  # array of Dec coordinates for plotting
	pointing_counter = 0    # the total number of pointings

	margin = 0.0           # determines edge-of-field margin, consistent with simdata

	# Everything in degrees:

	if ra_fov < 0.0:
		ra_fov = ra_fov * -1.0

	coords = generate_rectangle(ra_fov,dec_fov,centre_ra,centre_dec,theta)
	
	ra_fov = 1.3*np.abs((np.max(coords.T[0])-np.min(coords.T[0])))
	dec_fov = 1.3*np.abs((np.max(coords.T[1])-np.min(coords.T[1])))

	# Airy disk of a 12-m dish in degrees

	#pb_fwhm = 1.2 * 3.0e8/1.4e11/12.0*180.0/np.pi
	half_pb  = pb_fwhm * spacing

	ra_spacing = half_pb 
	dec_spacing = half_pb * 0.866025404  # cos 60 for hex

	n_rows = 1 + int(np.floor((dec_fov / dec_spacing) - 2.0 * margin / 0.866025404))

	float_cols = 1 + (ra_fov / ra_spacing) - (2.0 * margin)
	n_cols = int(np.floor(float_cols))

	if float_cols - n_cols >= 0.5 and n_rows > 1:
		even_cols = n_cols
		n_cols_min = 0.5 * (n_cols - 0.5)
	else:
		even_cols = n_cols -1
		n_cols_min = 0.5 * (n_cols - 1)

	current_dec = centre_dec + (0.5 * (n_rows-1) * dec_spacing)

	for i in range(0,n_rows):
		ra_spacing = half_pb / np.cos(current_dec*(np.pi/180.))

		if i % 2:
			ra_min = (centre_ra - (n_cols_min * ra_spacing))
			stop_col = n_cols
		else:
			ra_min = (centre_ra - ((n_cols_min - 0.5) * ra_spacing))
			stop_col = even_cols

		for j in range(0,stop_col):
			current_ra = ra_min + j * ra_spacing
			ra_mosaic_buffer.append(current_ra)
			dec_mosaic_buffer.append(current_dec)
			pointing_counter += 1
		current_dec = current_dec - dec_spacing
		
	temp_ra = []
	temp_dec =[]
	for i in range(len(ra_mosaic_buffer)):
		truth = determine_in_out(coords,[ra_mosaic_buffer[i],dec_mosaic_buffer[i]])
		if truth == True:
			temp_ra.append(ra_mosaic_buffer[i])
			temp_dec.append(dec_mosaic_buffer[i])
			
	ra_mosaic_buffer = temp_ra
	dec_mosaic_buffer = temp_dec
	# Write out a pointings file and generate a list of beams
	os.system('rm %s'%pointing_file_path)
	ptgfile = open(pointing_file_path,'w')

	print('#Epoch     RA   RA_hms       DEC  DEC_dms    RANGE',file=ptgfile)

	if ra_mosaic_buffer:
		for index in range(0, len(ra_mosaic_buffer)):
			c = SkyCoord(ra_mosaic_buffer[index],dec_mosaic_buffer[index],unit='deg')
			tmp_ra = str(ra_mosaic_buffer[index])
			ra_hms = str(int(c.ra.hms[0])).zfill(2)+'h'+str(int(c.ra.hms[1])).zfill(2)+"m"+str("%.12f"%(c.ra.hms[2])).zfill(15)+"s"
			tmp_dec = str(dec_mosaic_buffer[index])
			dec_dms = str(int(c.dec.dms[0])).zfill(2)+'d'+str(int(c.dec.dms[1])).zfill(2)+"m"+str("%.12f"%(c.dec.dms[2])).zfill(15)+"s"
			ptgstring = 'J2000 '+str(tmp_ra)+' '+ra_hms+' '+str(tmp_dec)+' '+dec_dms+' '+str(pb_fwhm)
			print(ptgstring, file=ptgfile)
	else:
		ptgstring = 'J2000 '+str(centre_ra)+' '+str(centre_dec)+' '+str(pb_fwhm)
		print(ptgstring,file=ptgfile)
	ptgfile.close()
	return coords
	#return ra_mosaic_buffer, dec_mosaic_buffer, half_pb, pointing_counter, coords

def hyp(co_a,co_b):
	return np.sqrt(np.abs(co_a[0]-co_b[0])**2. + np.abs(co_a[1]-co_b[1])**2)

def tri_area(co_a,co_b,co_c):
	##Get sizes of each side - assume co_d=[x_d,y_d]
	a = hyp(co_a,co_b)
	b = hyp(co_b,co_c)
	c = hyp(co_c,co_a)
	s = (a+b+c)/2
	return (s*(s-a)*(s-b)*(s-c)) ** 0.5

def rect_area(co_a,co_b,co_c):
	a = hyp(co_a,co_b)
	b = hyp(co_b,co_c)
	return a*b

def rotate(xy, theta):
	# https://en.wikipedia.org/wiki/Rotation_matrix#In_two_dimensions
	theta=theta*-1
	cos_theta, sin_theta = np.cos(theta), np.sin(theta)

	return np.array([
		xy[0] * cos_theta - xy[1] * sin_theta,
		xy[0] * sin_theta + xy[1] * cos_theta])


def translate(xy, offset):
	return np.array([xy[0] + offset[0], xy[1] + offset[1]])


def determine_in_out(coords,point):
	triangle_areas = []
	for i in range(len(coords)):
		if i == (len(coords)-1):
			triangle_areas.append(tri_area(coords[i],coords[0],point))
		else:
			triangle_areas.append(tri_area(coords[i],coords[i+1],point))
	tri_sum = np.sum(triangle_areas)
	rect_sum = rect_area(coords[0],coords[1],coords[2])
	return np.isclose(tri_sum,rect_sum)

def generate_central_wcs(crval, cdelt, crpix):
	# Create a new WCS object.  The number of axes must be set
	# from the start
	w = wcs.WCS(naxis=2)

	# Set up an "Airy's zenithal" projection
	# Vector properties may be set with Python lists, or Numpy arrays
	#CTYPE1  = projection
	#CRVAL1  = central position in degrees
	#CDELT1  = pixel demarcation
	#CRPIX1  = reference pixel
	#CUNIT1  = values of angle objects
	w.wcs.crpix = np.array(crpix).astype(int)
	w.wcs.cdelt = np.array(cdelt).astype(float)
	w.wcs.crval = np.array(crval).astype(float)
	w.wcs.ctype = ["RA---SIN", "DEC--SIN"]

	# Some pixel coordinates of interest.
	pixcrd = np.array([[-10, -10], [24, 38], [45, 98]], np.float_)

	# Convert pixel coordinates to world coordinates
	world = w.wcs_pix2world(pixcrd, 1)
	print(world)

	# Convert the same coordinates back to pixel coordinates.
	pixcrd2 = w.wcs_world2pix(world, 1)
	print(pixcrd2)

	# These should be the same as the original pixel coordinates, modulo
	# some floating-point error.
	assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6

	return w

def the_condition(xs,ys,condition):
	return xs.separation(ys).to(u.arcmin).value > condition ## arcmin separation to remove

def str_inp_convert(string):
	if ',' in str(string):
			ms2 = string.split(',')
			ms2_inp = ' '.join(ms2)
	return ms2_inp
def convert_frac_to_float(frac_str):
	try:
		return float(frac_str)
	except ValueError:
		num, denom = frac_str.split('/')
		try:
			leading, num = num.split(' ')
			whole = float(leading)
		except ValueError:
			whole = 0
		frac = float(num) / float(denom)
		return whole - frac if whole < 0 else whole + frac

# In[12]:

 
if str(inputs['custom_mosaic']) == '':
	cm=False
else:
	cm=True

if cm == False:
	pointing_centre = ast.literal_eval(inputs['field_centre'])
	pb_fwhm = np.load('%s/PB_fit.npy'%(inputs['output_path']))
	size = ast.literal_eval(inputs['mosaic_area']) 

	c = SkyCoord(pointing_centre[0], pointing_centre[1],unit=('hour','deg'))

	coords = mosaic_pointings_square(centre_ra=c.ra.deg,centre_dec=c.dec.deg,ra_fov=size[0]/np.cos(c.dec.rad), dec_fov=size[1],theta=0, pointing_file_path='%s/mosaic.csv'%(inputs['output_path']),pb_fwhm=pb_fwhm,spacing=float(inputs['mosaic_filling_factor']))

else:
	pb_fwhm = 1
	pointing_file_path='%s/mosaic.csv'%(inputs['output_path'])
	os.system('rm %s'%pointing_file_path)
	ptgfile = open(pointing_file_path,'w')
	with open('%s'%str(inputs['custom_mosaic'])) as f:
		lines = f.readlines()
	print('#Epoch     RA   RA_hms       DEC  DEC_dms    RANGE',file=ptgfile)
	mos_co = []
	for index in range(0, len(lines)):
			ra_mosaic_buffer = lines[index].strip('\n').split(' ')
			c = SkyCoord(ra_mosaic_buffer[0],ra_mosaic_buffer[1],unit=('hour','deg'))
			tmp_ra = str(c.ra.deg)
			ra_hms = str(int(c.ra.hms[0])).zfill(2)+'h'+str(int(c.ra.hms[1])).zfill(2)+"m"+str("%.12f"%(c.ra.hms[2])).zfill(15)+"s"
			tmp_dec = str(c.dec.deg)
			dec_dms = str(int(c.dec.dms[0])).zfill(2)+'d'+str(int(c.dec.dms[1])).zfill(2)+"m"+str("%.12f"%(c.dec.dms[2])).zfill(15)+"s"
			ptgstring = 'J2000 '+str(tmp_ra)+' '+ra_hms+' '+str(tmp_dec)+' '+dec_dms+' '+str(pb_fwhm)
			print(ptgstring, file=ptgfile)
			mos_co.append([float(tmp_ra),float(tmp_dec)])
	mos_co = np.array(mos_co)
	print(mos_co)
	c = SkyCoord(np.mean(mos_co[:,0]), np.mean(mos_co[:,1]),unit=('deg','deg'))
	ptgfile.close()

df= pd.read_csv('%s/mosaic.csv'%(inputs['output_path']), delim_whitespace=True)


w = generate_central_wcs([c.ra.deg,c.dec.deg],[1/60,1/60],[1,1])


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




