import os,sys
from collections import OrderedDict
import numpy as np
from scipy import constants as c
import pandas as pd
# CASA 6
import casatools
from casatasks import *
casalog.showconsole(True)
from astropy.io import fits
from astropy.nddata.utils import Cutout2D
from astropy.wcs import WCS
from astropy import units as u
from skimage.transform import resize, rotate, rescale
casa6=True

from astropy.utils.exceptions import AstropyWarning
import warnings
warnings.simplefilter('ignore', category=AstropyWarning)


def match_to_antenna_nos(evn_SEFD,msfile):
	tb = casatools.table()
	qa = casatools.quanta()
	me = casatools.measures()
	evn_SEFD_2 = {}
	evn_diams = {}
	tb.open('%s/ANTENNA'%msfile)
	x = tb.getcol('NAME')
	tb.close()
	print(evn_SEFD)
	for i,j in enumerate(x):
		evn_SEFD_2[i] = evn_SEFD[j][0]
		evn_diams[i] = evn_SEFD[j][1]
	return evn_SEFD_2, evn_diams

def P2R(radii, angles):
	return radii * np.exp(1j*angles)

def R2P(x):
	return np.abs(x), np.angle(x)

def calc_sefd(sefd1,sefd2,tint,dnu,sampeff):
	sefd_a = 1./np.sqrt(1./(float(sefd1)*float(sefd2)))
	#return(sefd_a)
	return (1./(float(sampeff)))*(sefd_a)*(1./(np.sqrt(tint*dnu)))

def add_noise(msfile,datacolumn,evn_SEFD,adjust_time=1.0):
	tb = casatools.table()
	qa = casatools.quanta()
	me = casatools.measures()
	tb.open('%s'%msfile, nomodify=True)
	data = tb.getcol('%s'%datacolumn)
	if datacolumn == 'CORRECTED_DATA':
		weightnames = 'WEIGHT'
	elif datacolumn == 'DATA':
		weightnames = 'SIGMA'
	else:
		raise TypeError
	weights = tb.getcol('%s'%weightnames)
	antenna1 = tb.getcol('ANTENNA1')
	antenna2 = tb.getcol('ANTENNA2')
	tint = np.average(tb.getcol('EXPOSURE'))
	tb.close()
	
	if adjust_time!=1.0:
		tint=tint*adjust_time

	tb.open('%s/SPECTRAL_WINDOW'%msfile, nomodify=True)
	chan_width=np.average(tb.getcol('CHAN_WIDTH'))
	print(chan_width,tint)
	tb.close()
	for i in range(len(antenna1)):
		sefd = calc_sefd(evn_SEFD[antenna1[i]],evn_SEFD[antenna2[i]],tint,chan_width,0.88)
		amps = np.random.normal(0.,sefd,np.shape(data[:,:,i]))
		phase = ((np.pi+np.pi)*np.random.random_sample(np.shape(data[:,:,i])))-np.pi
		data[:,:,i]=P2R(amps,phase)
		weights[:,i] = np.ones(weights[:,i].shape)/(sefd**2)

	tb.open('%s'%msfile, nomodify=False)
	tb.putcol('%s'%datacolumn,data)
	tb.putcol('%s'%weightnames,weights)
	tb.close()

def check_elevation(msfile,custom_xyz=False):
	tb = casatools.table()
	qa = casatools.quanta()
	me = casatools.measures()
	#measure = pm.measures()
	#tab = pt.table(msfile, readonly=True,ack=False)
	#field_tab = pt.table(tab.getkeyword('FIELD'),ack=False)
	#direction = np.squeeze(field_tab.getcol('PHASE_DIR'))
	
	tb.open(msfile, nomodify=True)
	time_unique = np.unique(tb.getcol('TIME'))
	tb.close()

	tb.open(msfile+'/FIELD',nomodify=True)
	direction = np.squeeze(tb.getcol('PHASE_DIR'))
	tb.close()

	tb.open(msfile+'/ANTENNA',nomodify=True)
	station_names = tb.getcol('NAME')
	if custom_xyz == True:
		if 'mosaic' in msfile:
			df = pd.read_csv('sims.itrf',delimiter=" ", header=None,names=['X', 'Y', 'Z', 'dish_diam', 'station', 'mount'],index_col=False)
		else:
			df = pd.read_csv('sims.itrf',delimiter=" ", header=None,names=['X', 'Y', 'Z', 'dish_diam', 'station', 'mount'],index_col=False)
		pos = np.vstack([df['X'].to_numpy(),df['Y'].to_numpy(),df['Z'].to_numpy()])
	else:
		pos = tb.getcol('POSITION')
	mount = tb.getcol('MOUNT')
	Nant = pos.shape[0]
	N = range(Nant)
	nbl = (Nant*(Nant-1))/2
	tb.close()

	ra = qa.quantity(direction[0], 'rad'); dec = qa.quantity(direction[1], 'rad')
	pointing = me.direction('j2000', ra, dec)
	start_time = me.epoch('utc', qa.quantity(time_unique[0], 's'))
	me.doframe(start_time)

	elevation_ant_matrix = np.zeros((Nant, time_unique.shape[0]))

	def antenna_elevation(antenna):
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
		return np.arcsin(np.sin(latitude)*np.sin(direction[1])+np.cos(latitude)*np.cos(direction[1]) *
						 np.cos(hour_angle))

	for i in range(Nant):
		elevation_ant_matrix[i] = antenna_elevation(i)
	return elevation_ant_matrix

def write_flag(msfile, elevation_limit, elevation, baseline_dict):
	""" flag data if below user-specified elevation limit """
	tb = casatools.table()
	qa = casatools.quanta()
	me = casatools.measures()
	tb.open(msfile, nomodify=True)
	flag = tb.getcol('FLAG').T
	tb.close()

	tb.open(msfile+'/ANTENNA',nomodify=True)
	station_names = tb.getcol('NAME')
	pos = tb.getcol('POSITION')
	mount = tb.getcol('MOUNT')
	Nant = pos.shape[0]
	tb.close()
	for a0 in range(Nant):
		for a1 in range(Nant):
			if a1 > a0:
				flag_mask = np.invert(((elevation[a1] > elevation_limit) & (elevation[a0] > elevation_limit)) > 0)
				#print(flag_mask.reshape((flag[:,:,baseline_dict[(a0, a1)]].shape, 1, 1)))
				print(flag_mask.reshape((flag_mask.shape[0], 1, 1)).shape)
				flag[baseline_dict[(a0, a1)]] = flag_mask.reshape((flag_mask.shape[0], 1, 1))
	if ("JB" in station_names) & ("M2" in station_names):
		flagdata(vis=msfile,mode='manual',antenna="JB&M2")
	tb.open(msfile, nomodify=False)
	tb.putcol("FLAG", flag.T)
	tb.close()

def make_baseline_dictionary(msfile):
	tb = casatools.table()
	qa = casatools.quanta()
	me = casatools.measures()
	tb.open(msfile, nomodify=True)
	A0 = tb.getcol('ANTENNA1')
	A1 = tb.getcol("ANTENNA2")
	ant_unique = np.unique(np.hstack((A0, A1)))
	tb.close()
	return dict([((x, y), np.where((A0 == x) & (A1 == y))[0])
				for x in ant_unique for y in ant_unique if y > x])


def calc_hpbw(x,diam,freq):
	x = (x/60.)*(np.pi/180.)
	wl = c.c/freq
	fwhm = wl/diam
	sigma = fwhm/(2*np.sqrt(2*np.log(2)))
	gaus = np.e**(-1*(x**2)/(2*sigma**2))
	return gaus
	
def calc_pb_corr(msfile,diam_ants,single_freq):
	tb = casatools.table()
	qa = casatools.quanta()
	me = casatools.measures()
	
	arcmin_off = float(msfile.split('_')[1]) + (float(msfile.split('_')[2].split('.ms')[0])/60.)

	tb.open('%s/SPECTRAL_WINDOW'%msfile)
	nspw = len(tb.getcol("MEAS_FREQ_REF"))
	chan_freqs = tb.getcol('CHAN_FREQ').T
	if single_freq!=False:
		freq=np.mean(chan_freqs)
		chan_freqs=[freq]*nspw
	diams_ants2 = {}
	for i in diams_ants.keys():
		pb_freq={}
		for j in range(nspw):
			pb_freq[str(j)] = np.sqrt(calc_hpbw(x=arcmin_off,diam=diam_ants[i],freq=chan_freqs[j]))
		print(pb_freq)
		diams_ants2[i] = pb_freq

	datacolumn='MODEL_DATA'
	tb.open('%s'%msfile, nomodify=True)
	data = tb.getcol('%s'%datacolumn)
	antenna1 = tb.getcol('ANTENNA1')
	antenna2 = tb.getcol('ANTENNA2')
	spw_id = tb.getcol('DATA_DESC_ID')
	tint = np.average(tb.getcol('EXPOSURE'))
	tb.close()
	for i in range(len(antenna1)):
		amps, phase = R2P(data[:,:,i])
		amps = amps*(diams_ants2[antenna1[i]][str(spw_id[i])]*diams_ants2[antenna2[i]][str(spw_id[i])])
		data[:,:,i]=P2R(amps,phase)
	tb.open('%s'%msfile, nomodify=False)
	tb.putcol('%s'%datacolumn,data)
	tb.close()
	return

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

def calc_pb_synthetic(msfile,diam_ants,single_freq,array):
	tb = casatools.table()
	qa = casatools.quanta()
	me = casatools.measures()
	degree_off = float(msfile.split('_')[1]) + (float(msfile.split('_')[2].split('.ms')[0])/60.)
	degree_off = degree_off/60.
	print('times')
	tb.open('%s'%msfile)
	t = tb.getcol('TIME')
	t = [np.min(t),np.max(t)]
	t = np.arange(t[0],t[1],3600)
	print(len(t))
	tb.close()
	print('parang')
	para_angle = parallacticAngle(msfile,t)
	print(np.max(para_angle),np.min(para_angle),np.shape(para_angle))
	if os.path.exists('../random_feed_rotation_%s.npy'%array):
		angle = np.load('../random_feed_rotation_%s.npy'%array)
		for i in range(np.shape(para_angle)[0]):
			para_angle[i,:] = (angle[i] + para_angle[i,:] + 180) % (2*180) - 180
	else:
		angle = np.random.default_rng().uniform(low=-180, high=180, size=(np.shape(para_angle)[0]))
		np.save('../random_feed_rotation_%s.npy'%array,angle)
		for i in range(np.shape(para_angle)[0]):
			para_angle[i,:] = (angle[i] + para_angle[i,:] + 180) % (2*180) - 180
	
	tb.open('%s/SPECTRAL_WINDOW'%msfile)
	nspw = len(tb.getcol("MEAS_FREQ_REF"))
	chan_freqs = tb.getcol('CHAN_FREQ').T
	tb.close()
	if single_freq!=False:
		freq=np.mean(chan_freqs)
		chan_freqs=[freq]*nspw

	datacolumn='MODEL_DATA'
	tb.open('%s'%msfile, nomodify=True)
	data = tb.getcol('%s'%datacolumn)
	times = tb.getcol('TIME')
	antenna1 = tb.getcol('ANTENNA1')
	antenna2 = tb.getcol('ANTENNA2')
	spw_id = tb.getcol('DATA_DESC_ID')
	tint = np.average(tb.getcol('EXPOSURE'))
	tb.close()
	print('rescale')
	for k in range(len(t)):
		diams_ants2 = {}
		for o,i in enumerate(diams_ants.keys()):
			pb_freq={}
			for j in range(nspw):
				print(k,o,j)
				pb_freq[str(j)] = rescale_synthetic_HPBW(header=[180.,60+degree_off],c_freq=chan_freqs[j],diameter=diam_ants[i],vmodel='%s_voltage_response_100.0m.fits'%msfile,a_term_upscale=8,phase_centre=[180.,60.],para_angle=para_angle[o,k])
				print('Ant %d (%s), time %.10e, freq %d, pbcor %.5f'%(o,i,t[k],chan_freqs[j],pb_freq[str(j)]))
				#pb_freq[str(j)] = np.sqrt(calc_hpbw(x=degree_off,diam=diam_ants[i],freq=chan_freqs[j]))
			diams_ants2[i] = pb_freq
	
		if k == (len(t)-1):
			sub_ant1 = antenna1[(times>=t[k])]
			sub_ant2 = antenna2[(times>=t[k])]
			sub_data = data[:,:,(times>=t[k])]
		else:
			sub_ant1 = antenna1[((times>=t[k])&(times<t[k+1]))]
			sub_ant2 = antenna2[((times>=t[k])&(times<t[k+1]))]
			sub_data = data[:,:,((times>=t[k])&(times<t[k+1]))]
		for i in range(len(sub_ant1)):
			amps, phase = R2P(sub_data[:,:,i])
			amps = amps*(diams_ants2[sub_ant1[i]][str(spw_id[i])]*diams_ants2[sub_ant2[i]][str(spw_id[i])])
			sub_data[:,:,i]=P2R(amps,phase)
		if k == (len(t)-1):
			data[:,:,(times>=t[k])] = sub_data
		else:
			data[:,:,((times>=t[k])&(times<t[k+1]))] = sub_data
	tb.open('%s'%msfile, nomodify=False)
	tb.putcol('%s'%datacolumn,data)
	tb.close()
	return

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

def add_pt_src(msfile,pt_flux):
	tb = casatools.table()
	qa = casatools.quanta()
	me = casatools.measures()
	cl = casatools.componentlist()
	tb.open(msfile+'/SOURCE')
	direc = tb.getcol('DIRECTION')
	direc = direc.T[0]
	tb.close()
	print('J2000 %srad %srad'%(direc[0],direc[1]))
	cl.addcomponent(flux=pt_flux, fluxunit='Jy',shape='point', dir='J2000 %srad %srad'%(direc[0],direc[1]))
	os.system('rm -r %s.cl'%msfile)
	cl.rename('%s.cl'%msfile)
	cl.close()
	ft(vis=msfile,complist='%s.cl'%msfile,usescratch=True)
	#uvsub(vis=msfile,reverse=True)
	os.system('rm -r %s.cl'%msfile)


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

ms = sys.argv[1]
imsize = int(sys.argv[2])
adjust_time = float(sys.argv[3])
band=str(sys.argv[4])

cell = str(sys.argv[5])
print('Clearing cal')
clearcal(vis=ms)
print('Write flag')
write_flag(ms,0,check_elevation(ms,custom_xyz=True),make_baseline_dictionary(ms))
print('Match antennae')
sefd_ants, diams_ants = match_to_antenna_nos(evn_SEFD[band],ms)
print('Add noise')
add_noise(msfile=ms,datacolumn='CORRECTED_DATA',evn_SEFD=sefd_ants,adjust_time=adjust_time)
print('Flag e-MERLIN')


os.system('rm -r %s_IM.*'%ms.split('.ms')[0])
tclean(vis=ms,
	   imagename='%s_IM'%ms.split('.ms')[0],
	   cell=cell,
	   imsize=[imsize,imsize],
	   deconvolver='clarkstokes',
	   niter=int(1e5),
	   nsigma=1.2,
	   usemask='user',
           pblimit=1e-10,
	   mask='circle[[%dpix, %dpix], 5pix]'%(imsize/2.,imsize/2.))
