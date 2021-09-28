import re
import sys
import traceback
import logging
import astropy.units as u
import numpy as np
from astropy.table import Table
from astropy import wcs
from astropy.coordinates import SkyCoord
def setup_logging_to_file(filename):
	logging.basicConfig( filename='./'+filename,
						 filemode='w',
						 level=logging.DEBUG,
						 format= '%(asctime)s - %(levelname)s - %(message)s',
					   )

def extract_function_name():
	"""Extracts failing function name from Traceback
	by Alex Martelli
	http://stackoverflow.com/questions/2380073/\
	how-to-identify-what-function-call-raise-an-exception-in-python
	"""
	tb = sys.exc_info()[-1]
	stk = traceback.extract_tb(tb, 1)
	fname = stk[0][3]
	return fname

def log_exception(e):
	logging.error(
	"Function {function_name} raised {exception_class} ({exception_docstring}): {exception_message}".format(
	function_name = extract_function_name(), #this is optional
	exception_class = e.__class__,
	exception_docstring = e.__doc__,
	exception_message = e.message))


def headless(inputfile):
	''' Parse the list of inputs given in the specified file. (Modified from evn_funcs.py)'''
	INPUTFILE = open(inputfile, "r")
	control = {}
	# a few useful regular expressions
	newline = re.compile(r'\n')
	space = re.compile(r'\s')
	char = re.compile(r'\w')
	comment = re.compile(r'#.*')
	# parse the input file assuming '=' is used to separate names from values
	for line in INPUTFILE:
		if char.match(line):
			line = comment.sub(r'', line)
			line = line.replace("'", '')
			(param, value) = line.split('=')
			param = newline.sub(r'', param)
			param = param.strip()
			param = space.sub(r'', param)
			value = newline.sub(r'', value)
			value = value.strip()
			valuelist = value.split(',')
			if len(valuelist) == 1:
				if valuelist[0] == '0' or valuelist[0]=='1' or valuelist[0]=='2':
					control[param] = int(valuelist[0])
				else:
					control[param] = str(valuelist[0])
			else:
				control[param] = ','.join(valuelist)
	return control

def generate_rectangle(x_fov,y_fov,RA_centre,Dec_centre,theta):
	fov = np.array([x_fov/2.,y_fov/2.])
	xy = np.array([fov,fov*np.array([-1,1]),fov*-1,fov*np.array([1,-1])]).T
	return translate(rotate(xy,theta*(np.pi/180.)),[RA_centre,Dec_centre]).T

def mosaic_pointings_square(centre_ra, centre_dec, centre_freq, ra_fov, dec_fov,theta, pointing_file_path):

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
	pb_fwhm=1./60.
	half_pb  = pb_fwhm * 0.6

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

	print('#Epoch     RA          DEC      RANGE',file=ptgfile)

	if ra_mosaic_buffer:
		for index in range(0, len(ra_mosaic_buffer)):
			tmp_ra = str(ra_mosaic_buffer[index])
			tmp_dec = str(dec_mosaic_buffer[index])
			ptgstring = 'J2000 '+str(tmp_ra)+' '+str(tmp_dec)+' '+str(pb_fwhm)
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

def filter_table(data,value):  ## Function to remove and filter the table of nearby matches
	result = []
	indices = []
	bad_result = []
	bad_indices = []
	repeats = []
	for i, element in enumerate(data):
		## Split the values into repeaters and non-repeaters
		if all(element.separation(other).to(u.arcmin).value > value for other in result):
			result.append(element)
			indices.append(i)
		else:
			bad_result.append(element)
			bad_indices.append(i)

	for j, index in enumerate(bad_indices):
		for k, goodind in enumerate(result):
			if bad_result[j].separation(result[k]).to(u.arcmin).value < value:
				if type(indices[k]) == int:
					indices[k] = [indices[k]]+[index]
				else:
					indices[k] = indices[k]+[index]
	return indices,result

def build_filtered_table(table, filter, filter_indices):
	if filter == 'False':
		return Table([table.ra.deg,table.dec.deg], names=('RA','DEC'))
	else: 
		df = Table([table.ra.deg,table.dec.deg], names=('RA','DEC'))
		RA = []
		DEC = []
		for i in range(len(filter_indices[0])):
			RA.append(np.average(df[filter_indices[0][i]]['RA']))
			DEC.append(np.average(df[filter_indices[0][i]]['DEC']))
		return Table([RA,DEC], names=('RA','DEC'))

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

def write_correlation_params(table,prefix):
	'''
	Function writes the correlation parameters in a $SOURCE vex format so that
	the correlator can read the phase centres in easily.
	e.g. def COSMOS-DEEP;
     source_name = COSMOS-DEEP;
	*    this source had calibrator code:  
	*    Center of COSMOS-deep
	     ra = 10h00m25.0000000s; dec =  02d33'00.000000"; ref_coord_frame = J2000;
	*    ra = 09h57m49.6816121s; dec =  02d47'25.904175"; ref_coord_frame = B1950;
	*    ra = 10h01m24.8107239s; dec =  02d27'22.307892"; ref_coord_frame = Date;
	enddef;
	'''
	correlation_string = ['$SOURCE;']
	for i in range(len(table['RA'])):
		sig_fig = len(str(len(table['RA'])))
		c = SkyCoord(table['RA'][i],table['DEC'][i],unit=('deg','deg'))
		sky_string = c.to_string('hmsdms').split(' ')
		RA_string = sky_string[0]
		Dec_string = sky_string[1].replace('m',"'").replace('s','\"')
		if Dec_string[0] == '+':
			Dec_string = Dec_string.replace('+',' ')
		correlation_string.append('def %s%s;'%(prefix[0:(8-sig_fig)],'{0:0{width}}'.format(i, width=sig_fig)))
		correlation_string.append('    source_name = %s%s;'%(prefix[0:(8-sig_fig)],'{0:0{width}}'.format(i, width=sig_fig)))
		correlation_string.append('    ra = %s; dec = %s; ref_coord_frame = J2000;' % (RA_string,Dec_string))
		correlation_string.append('enddef;')
	return correlation_string
