import re, os, inspect, sys
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
import matplotlib.pyplot as plt
from astropy import wcs
import astropy.units as u
from astropy.io import fits
from matplotlib.colors import SymLogNorm
import pandas as pd

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

def write_hpc_headers(step,params):
	func_name = inspect.stack()[0][3]
	print(params)
	hpc_opts = {}
	hpc_opts['job_manager'] = params['job_manager']
	hpc_opts['job_name'] = 'vlbisched_%s'%step
	hpc_opts['email_progress'] = params["email_progress"] 
	hpc_opts['hpc_account'] = params['HPC_project_code']
	hpc_opts['error'] = step

	if ((hpc_opts['job_manager'] == 'pbs')|(hpc_opts['job_manager'] == 'bash')|(hpc_opts['job_manager'] == 'slurm')):
		pass
	#else:
		#sys.exit()

	for i in ['partition','walltime','nodetype','nodes','cpus','mpiprocs']:
		hpc_opts[i] = params['%s'%i]


	hpc_dict = {'slurm':{
					 'partition'     :'#SBATCH --partition=%s'%hpc_opts['partition'],
					 'nodetype'      :'',
					 'cpus'          :'#SBATCH --tasks-per-node %s'%hpc_opts['cpus'], 
					 'nodes'         :'#SBATCH -N %s-%s'%(hpc_opts['nodes'],hpc_opts['nodes']),
					 'mpiprocs'      :'', 
					 'walltime'      :'#SBATCH --time=%s'%hpc_opts['walltime'],
					 'job_name'      :'#SBATCH -J %s'%hpc_opts['job_name'],
					 'hpc_account'   :'#SBATCH --account %s'%hpc_opts['hpc_account'],
					 'email_progress':'#SBATCH --mail-type=BEGIN,END,FAIL\n#SBATCH --mail-user=%s'%hpc_opts['email_progress'],
					 'error':'#SBATCH -o logs/%s.sh.stdout.log\n#SBATCH -e logs/%s.sh.stderr.log'%(hpc_opts['error'],hpc_opts['error'])
					},
				'pbs':{
					 'partition'     :'#PBS -q %s'%hpc_opts['partition'],
					 'nodetype'      :'',
					 'cpus'          :'#PBS -l select=%s:ncpus=%s:mpiprocs=%s:nodetype=%s'%(hpc_opts['nodes'],hpc_opts['cpus'],hpc_opts['mpiprocs'],hpc_opts['nodetype']), 
					 'nodes'         :'',
					 'mpiprocs'      :'', 
					 'walltime'      :'#PBS -l walltime=%s'%hpc_opts['walltime'],
					 'job_name'      :'#PBS -N %s'%hpc_opts['job_name'],
					 'hpc_account'   :'#PBS -P %s'%hpc_opts['hpc_account'],
					 'email_progress':'#PBS -m abe -M %s'%hpc_opts['email_progress'],
					 'error':'#PBS -o logs/%s.sh.stdout.log\n#PBS -e logs/%s.sh.stderr.log'%(hpc_opts['error'],hpc_opts['error'])
					},
				'bash':{
					 'partition'     :'',
					 'nodetype'      :'',
					 'cpus'          :'', 
					 'nodes'         :'',
					 'mpiprocs'      :'', 
					 'walltime'      :'',
					 'job_name'      :'',
					 'hpc_account'   :'',
					 'email_progress':'',
					 'error':''
					}
				}

	hpc_header= ['#!/bin/bash']

	if step == 'run_mosaic_sims':
		file = open("%s/target_files.txt"%params['global']['cwd'], "r")
		nonempty_lines = [line.strip("\n") for line in file if line != "\n"]
		line_count = len(nonempty_lines)
		file.close()
		if params[step]['hpc_options']['max_jobs'] == -1:
			tasks = '0-'+str(line_count-1)
		else:
			if (line_count-1) > params[step]['hpc_options']['max_jobs']:
				tasks = '0-'+str(line_count-1)+'%'+str(params[step]['hpc_options']['max_jobs'])
			else:
				tasks = '0-'+str(line_count-1)
		hpc_dict['slurm']['array_job'] = '#SBATCH --array='+tasks
		hpc_dict['pbs']['array_job'] = '#PBS -t '+tasks
		hpc_dict['bash']['array_job'] = ''
		hpc_opts['array_job'] = -1

	hpc_job = hpc_opts['job_manager']
	for i in hpc_opts.keys():
		if i != 'job_manager':
			if hpc_opts[i] != '':
				if hpc_dict[hpc_opts['job_manager']][i] !='':
					hpc_header.append(hpc_dict[hpc_job][i])


	with open('job_%s.%s'%(step,hpc_job), 'w') as filehandle:
		for listitem in hpc_header:
			filehandle.write('%s\n' % listitem)