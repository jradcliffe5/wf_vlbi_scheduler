from simms import simms
import os, glob, sys
import numpy as np
from datetime import datetime, timedelta

#INPUT  = "input"
#OUTPUT = "output"
#MSDIR  = "msdir"
#singularity_image_dir=os.environ["STIMELA_IMAGES_PATH"]
#stimela.register_globals()
#SKYMODEL   = "blank_sky.txt"
import sys
arrays = sys.argv[2:]
arrays = "_".join(arrays)

if sys.argv[1] == 'single':
	os.system('rm -r %s.ms'%arrays)
	MS='%s.ms'%arrays
	simms.create_empty_ms(
	msname=MS,
	label=None,
	tel="EVN",
	pos="%s_vlapos_sims.itrf"%arrays,
	pos_type='ascii',
	ra="12h00m00s",
	dec="60d00m00s",
	synthesis=12,
	scan_length=[1],
	dtime=10,
	freq0="1.536GHz",#1536000000.0,
	dfreq="4MHz",
	nchan="32",
	stokes='RR RL LR LL',
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
	os.system('rm -r mosaic_ms')
	os.system('mkdir mosaic_ms')
	with open('mosaic.csv') as f:
		lines = f.readlines()
	direction = []
	for i in lines:
		if not i.startswith('#'):
			line = i.split(" ")
			direction.append([line[2],line[4]])
	total_time = 200.0
	for i in range(len(direction)):
		print(direction[i])
		dt = datetime.strptime('20 Sep 2021', '%d %b %Y') #+ timedelta(hours=2/60+(0*total_time/float(len(direction))))
		MS='mosaic_ms/%s_mosaic_%s.ms'%(arrays,i)
		simms.create_empty_ms(
		msname=MS,
		label=None,
		tel="EVN",
		pos="%s_vlapos_sims.itrf"%arrays,
		pos_type='ascii',
		ra=direction[i][0],
		dec=direction[i][1],
		direction=None,
		synthesis=total_time/float(len(direction)),
		scan_length=[1],
		dtime=10,
		freq0="1.536GHz",#1536000000.0,
		dfreq="4MHz",
		nchan="32",
		stokes='RR RL LR LL',
		setlimits=False,
		elevation_limit=0,
		shadow_limit=0,
		outdir="./",
		nolog=False,
		coords='itrf',
		lon_lat=None,
		noup=False,
		nbands=1,
		date="UTC,%s"%dt.strftime('%Y/%m/%d/%H:%M:%S'),
		fromknown=False,
		feed='perfect R L',
		scan_lag=0,
		auto_corr=False,
		optimise_start=None
		)
else:
	print('Incorrect entry ... exiting')
	sys.exit()
