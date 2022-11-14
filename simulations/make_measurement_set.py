from simms import simms
import os, glob, sys, ast
import numpy as np
from datetime import datetime, timedelta
from simulator_functions import headless

#INPUT  = "input"
#OUTPUT = "output"
#MSDIR  = "msdir"
#singularity_image_dir=os.environ["STIMELA_IMAGES_PATH"]
#stimela.register_globals()
#SKYMODEL   = "blank_sky.txt"
import sys
inputs = headless(sys.argv[2])
output = str(inputs['output_path'])

## Calculate datarates
data_rate = float(inputs['data_rate'])
npols = float(inputs['npols'])
if npols >=2.:
	bit = 2.
	if npols == 4.:
		stokes ='RR RL LR LL'
	else:
		stokes='RR LL'
else:
	bit = 1.
	stokes = 'RR'
try: 
	bw = float(inputs['bandwidth'])
except:
	bw = data_rate/bit/4.
freq0 = float(inputs['obs_freq'])-(bw/2000.)

if inputs['mosaic'] == "True":
	tos = 12
else:
	tos = float(inputs['total_time_on_source'])

if str(inputs['wide_field_ITRF']) == 'True':
	itrf="%s/vlapos_sims.itrf"%output
else:
	itrf="%s/sims.itrf"%output

pointing_centre = ast.literal_eval(inputs['field_centre'])

if sys.argv[1] == 'single':
	os.system('rm -r %s/single_pointing.ms'%output)
	tos = np.round(tos,5)
	MS='%s/single_pointing.ms'%output
	simms.create_empty_ms(
	msname=MS,
	label=None,
	tel="EVN",
	pos=itrf,
	pos_type='ascii',
	ra=pointing_centre[0],
	dec=pointing_centre[1],
	synthesis=tos,
	scan_length=[tos],
	dtime=30,
	freq0="%.8fGHz"%freq0,#1536000000.0,
	dfreq="%.8fMHz"%(bw/64.),
	nchan="64",
	stokes=stokes,
	setlimits=False,
	elevation_limit=0,
	shadow_limit=0,
	outdir="./",
	nolog=False,
	coords='itrf',
	lon_lat=None,
	noup=False,
	nbands=1,
	direction=[],
	date=None,
	fromknown=False,
	feed='perfect R L',
	scan_lag=0,
	auto_corr=False,
	optimise_start=None
	)
elif sys.argv[1] == 'mosaic':
	with open('mosaic.csv') as f:
		lines = f.readlines()
	direction = []
	for i in lines:
		if not i.startswith('#'):
			line = i.split(" ")
			direction.append([line[2],line[4]])
	tos = float(inputs['total_time_on_source'])
	total_time = tos
	print(total_time, total_time/float(len(direction)) , len(direction))
	synthesis = np.round(total_time/float(len(direction)),3)
	for i in range(len(direction)):
		os.system('rm -r %s/mosaic_%s.ms'%(output,i))
		dt = datetime.strptime('20 Sep 2021', '%d %b %Y') #+ timedelta(hours=2/60+(0*total_time/float(len(direction))))
		MS='%s/mosaic_%s.ms'%(output,i)
		simms.create_empty_ms(
		msname=MS,
		label=None,
		tel="EVN",
		pos=itrf,
		pos_type='ascii',
		ra=direction[i][0],
		dec=direction[i][1],
		direction=None,
		synthesis=synthesis,
		scan_length=[synthesis],
		dtime=10,
		freq0="%.8fGHz"%freq0,#1536000000.0,
		dfreq="%.8fMHz"%(bw/32.),
		nchan="32",
		stokes=stokes,
		setlimits=False,
		elevation_limit=0,
		shadow_limit=0,
		outdir="./",
		nolog=False,
		coords='itrf',
		lon_lat=None,
		noup=False,
		nbands=1,
		date=None,
		fromknown=False,
		feed='perfect R L',
		scan_lag=0,
		auto_corr=False,
		optimise_start=None
		)
else:
	print('Incorrect entry ... exiting')
	sys.exit()
