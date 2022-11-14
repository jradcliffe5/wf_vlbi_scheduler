import os, sys
from collections import OrderedDict
import numpy as np
from scipy import constants as c
import pickle

try:
	# CASA 6
	import casatools
	from casatasks import *
	casalog.showconsole(True)
	from astropy.io import fits
	from astropy.nddata.utils import Cutout2D
	from astropy.wcs import WCS
	from astropy import units as u
	from astropy.coordinates import SkyCoord
	from skimage.transform import resize, rotate, rescale
	from astropy.utils.exceptions import AstropyWarning
	import warnings
	from scipy import constants as c
	warnings.simplefilter('ignore', category=AstropyWarning)
	casa6=True
except:
	# CASA 5
	from casac import casac as casatools
	from taskinit import casalog
	import pyfits as fits
	casa6=False

try:
	# Python 2
	from StringIO import StringIO
except:
	# Python 3
	from io import StringIO

def parallacticAngle(msfile,times):
	#measure = pm.measures()
	#tab = pt.table(msfile, readonly=True,ack=False)
	#field_tab = pt.table(tab.getkeyword('FIELD'),ack=False)
	#direction = np.squeeze(field_tab.getcol('PHASE_DIR'))
	tb = casatools.table()
	qa = casatools.quanta()
	me = casatools.measures()
	time_unique = times

	tb.open(msfile+'/FIELD',nomodify=True)
	direction = np.squeeze(tb.getcol('PHASE_DIR'))
	tb.close()

	tb.open(msfile+'/ANTENNA',nomodify=True)
	station_names = tb.getcol('NAME')
	pos = tb.getcol('POSITION').T
	mount = tb.getcol('MOUNT')
	Nant = pos.shape[0]
	N = range(Nant)
	nbl = (Nant*(Nant-1))/2
	tb.close()

	ra = qa.quantity(direction[0], 'rad'); dec = qa.quantity(direction[1], 'rad')
	pointing = me.direction('j2000', ra, dec)
	start_time = me.epoch('utc', qa.quantity(time_unique[0], 's'))
	me.doframe(start_time)

	parallactic_ant_matrix = np.zeros((Nant, time_unique.shape[0]))

	def antenna_para(antenna):
		x = qa.quantity(pos[antenna, 0], 'm')
		y = qa.quantity(pos[antenna, 1], 'm')
		z = qa.quantity(pos[antenna, 2], 'm')
		position = me.position('wgs84', x, y, z)
		me.doframe(position)
		sec2rad = 2 * np.pi / (24 * 3600.)
		hour_angle = me.measure(pointing, 'HADEC')['m0']['value'] +\
					 (time_unique-time_unique.min()) * sec2rad
		earth_radius = 6371000.0
		latitude = np.arcsin(pos[antenna, 2]/earth_radius)
		return np.arctan2(np.sin(hour_angle)*np.cos(latitude), (np.cos(direction[1])*np.sin(latitude)-np.cos(hour_angle)*np.cos(latitude)*np.sin(direction[1])))

	for i in range(Nant):
		if mount[i] == 'EQUATORIAL':
			parallactic_ant_matrix[i] = np.zeros(time_unique.shape)
		else:
			parallactic_ant_matrix[i] = antenna_para(i)*(180./np.pi)
	return parallactic_ant_matrix

def match_to_antenna_nos(evn_SEFD,msfile):
	evn_SEFD_2 = {}
	evn_diams = {}
	tb.open('%s/ANTENNA'%msfile)
	x = tb.getcol('NAME')
	tb.close()
	for i,j in enumerate(x):
		evn_SEFD_2[i] = evn_SEFD[j][0]
		evn_diams[i] = evn_SEFD[j][1]
	return evn_SEFD_2, evn_diams

def calc_hpbw(x,diam,freq):
	x = (x/60.)*(np.pi/180.)
	wl = c.c/freq
	fwhm = wl/diam
	sigma = fwhm/(2*np.sqrt(2*np.log(2)))
	gaus = np.e**(-1*(x**2)/(2*sigma**2))
	return gaus

def same_dist_elems(arr):
	diff = arr[1] - arr[0]
	for x in range(1, len(arr) - 1):
		if arr[x + 1] - arr[x] != diff:
			return False
	return True

def rescale_synthetic_HPBW(header,c_freq,diameter,vmodel,a_term_upscale,phase_centre,para_angle):
	tb = casatools.table()
	qa = casatools.quanta()
	me = casatools.measures()
	try:
		### Set sizes
		size = np.array([header['NAXIS1'],header['NAXIS2']])
		cdelt = np.array([header['CDELT1']*a_term_upscale,header['CDELT2']*a_term_upscale])
		centre = np.array([header['CRVAL1'],header['CRVAL2']])
		total_size = size*cdelt
		#print('Header found')
		ret_singular = False
	except:
		#print('No header found reverting to some default values')
		size = [512,512]
		cdelt = [1.388888888889e-07,1.388888888889e-07]
		centre = [header[0],header[1]]
		ret_singular = True
	
	### Open voltage model
	#print('open fits')
	hdu = fits.open(vmodel)
	vhead = hdu[0].header
	#print('reorder fits')
	vdata = hdu[0].data.squeeze()
	vdata = vdata.byteswap().newbyteorder()
	hdu.close()
	
	### Rotate data
	#try:
	#print('sklearn rotate')
	ct = np.where(vdata==vdata.max())
	#print('rotate')
	#print(vdata.shape,para_angle,ct[0][0],ct[1][0])
	vdata = rotate(image=vdata,angle=para_angle,center=[ct[0][0],ct[1][0]])
	#except:
	#from scipy.ndimage import rotate
	#print('scipy rotate')
	#vdata = rotate(input=vdata,angle=para_angle,reshape=False)
	
	### Rescale model to diameter
	vhead['CDELT1'] = vhead['CDELT1']*(vhead['DIAMETER']/diameter)
	vhead['CDELT2'] = vhead['CDELT2']*(vhead['DIAMETER']/diameter)
	
	### Rescale model to frequency
	vhead['CDELT1'] = vhead['CDELT1']*(vhead['CRVAL3']/c_freq)
	vhead['CDELT2'] = vhead['CDELT2']*(vhead['CRVAL3']/c_freq)
	
	### Set phase centres
	vhead['CRVAL1'] = phase_centre[0]
	vhead['CRVAL2'] = phase_centre[1]
	
	#print('wcs')
	### Set cutout re-sampling
	wcs = WCS(vhead,naxis=2)
	#print('wcs convert')
	centre_p = wcs.all_world2pix(centre[0]*u.deg,centre[1]*u.deg,0)
	
	scale = vhead['CDELT1']/np.abs(cdelt[0])
	
	## Add hard total pixel scaling factor 
	# (bigger machines can have > 2e4 but probably overkill)
	sz=2000
	while scale*sz > 2e4:
		sz = int(sz-2)
		if sz < 2:
			print("error size less than 2")
			sys.exit()
	hdu_cut = Cutout2D(data=vdata,wcs=wcs,\
			 position=[centre_p[0],centre_p[1]],\
			 size=[sz,sz])
	
	
	## Adjust scales to take this into account
	vhead = hdu_cut.wcs.to_header()
	vhead['CDELT1'] = vhead['CDELT1']/scale
	vhead['CDELT2'] = vhead['CDELT2']/scale
	vhead['CRPIX1'] = vhead['CRPIX1']*scale
	vhead['CRPIX2'] = vhead['CRPIX2']*scale
	
	data = rescale(hdu_cut.data,scale)
	wcs = WCS(vhead,naxis=2)
	centre_p = wcs.all_world2pix(centre[0]*u.deg,centre[1]*u.deg,0)
	hdu_cut = Cutout2D(data=data,
						wcs=wcs,\
						position=[centre_p[0],centre_p[1]],\
						 size=size)

	if ret_singular == False:
		return hdu_cut.data
	else:
		hc = hdu_cut.data.shape
		return hdu_cut.data[int(hc[0]/2),int(hc[1]/2)]



evn_SEFD = {'L':{
				'Ef':[19,76],
				'Tm65':[39,65],
				'Jb1':[65,67],
				'W1':[560,25],
				'On':[350,25],
				'Mc':[700,32],
				'Tr':[300,32],
				'Nt':[740,25],
				'Sv':[360,32],
				'Bd':[330,32],
				'Zc':[300,32],
				'Ur':[300,25],
				'Cm':[212,32],
				'Da':[450,25],
				'Kn':[400,25],
				'Pi':[450,25],
				'De':[350,25],
				'Sh':[670,25],
				'Ir':[3600,25],
				'Jb2':[320,25],
				'Ys':[160,25],
				'Sc':[365,25],
				'Hn':[365,25],
				'Nl':[365,25],
				'Fd':[365,25],
				'La':[365,25],
				'Kp':[365,25],
				'Pt':[365,25],
				'Ov':[365,25],
				'Br':[365,25],
				'Mk':[365,25]},
			'C':{'Ef':[20,76],
				'Tm65':[39,65],
				'Jb1':[40,67],
				'W1':[840,25],
				'On':[600,25],
				'Mc':[170,32],
				'Tr':[220,32],
				'Nt':[260,25],
				'Sv':[250,32],
				'Bd':[200,32],
				'Zc':[400,32],
				'Ur':[200,25],
				'Cm':[136,32],
				'Da':[325,25],
				'Kn':[325,25],
				'Pi':[325,25],
				'De':[1000,25],
				'Sh':[720,25],
				'Ir':[430,25],
				'Jb2':[320,25],
				'Ys':[160,25]}
			}

if int(sys.argv[3]) == 1:
	if os.path.exists('../D_eff_errs.pkl'):
		infile = open('../D_eff_errs.pkl','rb')
		evn_SEFD = pickle.load(infile)
	else:
		outfile = open('../D_eff_errs.pkl','wb')
		for i in evn_SEFD.keys():
			pc_err = 10.
			evn_SEFD[i][1] = np.random.normal(loc=evn_SEFD[i][1], scale=(evn_SEFD[i][1]*(pc_err/100.)), size=1)
		pickle.dump(evn_SEFD,outfile)
		outfile.close()

ms = sys.argv[1]
band = str(sys.argv[5])
if int(sys.argv[2]) == 1:
	do_all_freqs = True
else:
	do_all_freqs = False
if int(sys.argv[4]) >= 1:
	do_time_rotation = True
	if int(sys.argv[4]) == 2:
		dont_rotate = True
	else:
		dont_rotate = False
else:
	do_time_rotation = False
	dont_rotate = True

tb = casatools.table()
qa = casatools.quanta()
me = casatools.measures()
try:
	ang_off = float(ms.split('_')[1]) + float(ms.split('_')[2].split('.ms')[0])/60.
	degree_off = ang_off/60.
	diffcorr = False
except:
	diffcorr = True
evn_sefd, evn_diams = match_to_antenna_nos(evn_SEFD[band],ms)
nants = len(evn_diams.values())

os.system('rm -r %s.mask'%ms)
ia = casatools.image()
ia.open('%s_IM.image'%ms.split('.ms')[0])
imsize=ia.shape()[0] 
ia.close()
makemask(inpimage='%s_IM.image'%ms.split('.ms')[0],mode='copy',inpmask='circle[[%dpix, %dpix], 5pix]'%(imsize/2.,imsize/2.),output=ms+'.mask')
exportfits(imagename='%s_IM.image'%ms.split('.ms')[0],fitsimage='%s_IM.image.fits'%ms.split('.ms')[0],overwrite=True)
hdu1 = fits.open('%s_IM.image.fits'%ms.split('.ms')[0])
header=hdu1[0].header

tb.open('%s/SPECTRAL_WINDOW'%ms)
nspw = len(tb.getcol("MEAS_FREQ_REF"))
chan_freq = tb.getcol('CHAN_FREQ').T.flatten()
tb.close()

tb.open('%s/FIELD'%ms)
field_dir = tb.getcol('PHASE_DIR').squeeze()

if same_dist_elems(chan_freq) == True:
	freq0 = chan_freq[0]
	chan_width = chan_freq[1]-chan_freq[0]
	bw = chan_freq[-1]-chan_freq[0] + np.abs(chan_width)
	freqs = np.linspace(freq0,freq0+((len(chan_freq)-1)*chan_width),len(chan_freq))
else:
	sys.exit()

header['CDELT2'] = 3*header['CDELT2']
header['CDELT1'] = 3*header['CDELT1']
header['CTYPE3']  = 'MATRIX'                                                            
header['CRPIX3']  =                   1.
header['CRVAL3']  =                   0.   
header['CTYPE4']  = 'ANTENNA'                                                            
header['CRPIX4']  =                   1.
header['CRVAL4']  =                   0. 
header['CTYPE5']  = 'FREQ' 
header['CRPIX5']  =   1. 
if do_all_freqs == True:                                                            
	header['CRVAL5']  =   freq0
	header['CDELT5']  =   chan_width
else:                                                           
	header['CRVAL5']  =   freq0                                                  
	header['CDELT5']  =   bw
header['CUNIT5']  = 'Hz      '
header['CTYPE6']  = 'TIME'
header['CRPIX6']  =                   1.
if do_time_rotation == True: 
	header['CRVAL6']  =    t_range[0]+interval/2. 
	header['CDELT6']  =        float(interval)
else:
	header['CRVAL6']  =                   1.
	header['CDELT6']  =             1.0E+20
print(header['CRVAL5'])

if do_time_rotation == True:
	print('rotating')
	interval=3600
	tb.open('%s'%ms)
	t = tb.getcol('TIME')
	t_range = [np.min(t),np.max(t)]
	t = np.arange(t_range[0],t_range[1],interval)
	tb.close()
	para_angle= parallacticAngle(ms,t)
	if dont_rotate == True:
		para_angle = para_angle*0
	if os.path.exists('../random_feed_rotation_%s.npy'%ms.split('_')[0]):
		angle = np.load('../random_feed_rotation_%s.npy'%ms.split('_')[0])
		for i in range(np.shape(para_angle)[0]):
			para_angle[i,:] = (angle[i] + para_angle[i,:] + 180) % (2*180) - 180
	else:
		sys.exit()

if do_all_freqs == True:
	if do_time_rotation == True:
		x = np.zeros([len(t),len(freqs),nants,4,header['NAXIS1'],header['NAXIS2']],dtype=np.float32)
		if dont_rotate == True:
			for i in range(nants):
				for k,j in enumerate(freqs):
					corr = rescale_synthetic_HPBW(header=[180.,60+degree_off],c_freq=j,diameter=evn_diams[i],vmodel='%s_voltage_response_100.0m.fits'%ms,a_term_upscale=8,phase_centre=[180.,60.], para_angle=para_angle[i,0])
					x[:,k,i,0,:,:]= corr
					x[:,k,i,2,:,:]= corr
					print('Ant %d, freq %d, pbcor %.5f'%(i,j,corr))
		else:
			for o in range(len(t)):
				for i in range(nants):
					for k,j in enumerate(freqs):
						corr = rescale_synthetic_HPBW(header=[180.,60+degree_off],c_freq=j,diameter=evn_diams[i],vmodel='%s_voltage_response_100.0m.fits'%ms,a_term_upscale=8,phase_centre=[180.,60.], para_angle=para_angle[i,o])
						x[o,k,i,0,:,:]= corr
						x[o,k,i,2,:,:]= corr
						print('Ant %d, time %.10e, freq %d, pbcor %.5f'%(i,t[o],j,corr))
	else:
		x = np.zeros([1,len(freqs),nants,4,header['NAXIS1'],header['NAXIS2']],dtype=np.float32)
		if diffcorr:
			w = WCS(header,naxis=2)
			phasecentre = SkyCoord(field_dir[0],field_dir[1],unit=('rad','rad'))
			xx = np.arange(1,int(header['NAXIS1'])+1,1)
			yy = np.arange(1,int(header['NAXIS2'])+1,1)
			xv, yv = np.meshgrid(xx, yy, sparse=True, indexing='xy')
			c0 = w.celestial.wcs_pix2world(xv,yv,1)
			c1 = SkyCoord(c0[0],c0[1],unit=('deg','deg'))
			pc = phasecentre
			ang_off = np.float64(c1.separation(pc).rad)
			wl = c.c/(header['CRVAL5']+(header['CDELT5']/2.))
		for i in range(nants):
			for k,j in enumerate(freqs):
				corr = np.sqrt(calc_hpbw(ang_off,evn_diams[i],j))
				x[:,k,i,0,:,:] = corr
				x[:,k,i,2,:,:] = corr
				print('Ant %d, freq %d, pbcor %.5f'%(i,j,corr))
else:
	c_freq = np.mean(chan_freq.flatten())
	if do_time_rotation == True:
		x = np.zeros([len(t),1,nants,4,header['NAXIS1'],header['NAXIS2']],dtype=np.float32)
		if dont_rotate == True:
			for i in range(nants):
				corr = rescale_synthetic_HPBW(header=[180.,60+degree_off],c_freq=c_freq,diameter=evn_diams[i],vmodel='%s_voltage_response_100.0m.fits'%ms,a_term_upscale=8,phase_centre=[180.,60.], para_angle=para_angle[i,0])
				x[:,:,i,0,:,:]= corr
				x[:,:,i,2,:,:]= corr
				print('Ant %d, freq %d, pbcor %.5f'%(i,c_freq,corr))
		else:
			for o in range(len(t)):
				for i in range(nants):
					corr = rescale_synthetic_HPBW(header=[180.,60+degree_off],c_freq=c_freq,diameter=evn_diams[i],vmodel='%s_voltage_response_100.0m.fits'%ms,a_term_upscale=8,phase_centre=[180.,60.], para_angle=para_angle[i,o])
					x[o,:,i,0,:,:]= corr
					x[o,:,i,2,:,:]= corr
					print('Ant %d, time %.10e, freq %d, pbcor %.5f'%(i,t[o],c_freq,corr))

	else:
		x = np.zeros([1,1,nants,4,header['NAXIS1'],header['NAXIS2']],dtype=np.float32)
		if diffcorr:
			w = WCS(header,naxis=2)
			phasecentre = SkyCoord(field_dir[0],field_dir[1],unit=('rad','rad'))
			xx = np.arange(1,int(header['NAXIS1'])+1,1)
			yy = np.arange(1,int(header['NAXIS2'])+1,1)
			xv, yv = np.meshgrid(xx, yy, sparse=True, indexing='xy')
			c0 = w.celestial.wcs_pix2world(xv,yv,1)
			c1 = SkyCoord(c0[0],c0[1],unit=('deg','deg'))
			pc = phasecentre
			ang_off = np.float64(c1.separation(pc).arcmin)
			wl = c.c/(header['CRVAL5']+(header['CDELT5']/2.))
		for i in range(nants):
			corr = np.sqrt(calc_hpbw(ang_off,evn_diams[i],c_freq))
			x[:,:,i,0,:,:] = corr
			x[:,:,i,2,:,:] = corr



hdu = fits.PrimaryHDU(x,header)
header = hdu.header

if dont_rotate == True:
	hdu.writeto('%s_pb_flat_norotate.fits'%(ms),overwrite=True)
else:
	hdu.writeto('%s_pb_flat.fits'%(ms),overwrite=True)

if dont_rotate == True:
	preamble = ['# a term corrections', '# This is a test parset, comments are possible like this',
			'aterms = [diagonal]', 'diagonal.images = [ %s ]'%'%s_pb_flat_norotate.fits'%(ms), 'diagonal.tukeywindow = false']
	with open('%s_aterm_norotate_config.txt'%(ms), 'w') as f:
		for item in preamble:
			f.write("%s\n" % item)
else:
	preamble = ['# a term corrections', '# This is a test parset, comments are possible like this',
			'aterms = [diagonal]', 'diagonal.images = [ %s ]'%'%s_pb_flat.fits'%(ms), 'diagonal.tukeywindow = false']
	with open('%s_aterm_config.txt'%(ms), 'w') as f:
		for item in preamble:
			f.write("%s\n" % item)


