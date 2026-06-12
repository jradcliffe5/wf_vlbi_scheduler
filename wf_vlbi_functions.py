import re, os
import sys
import json
import traceback
import logging
import astropy.units as u
import numpy as np
from astropy.table import Table
from astropy import wcs
from astropy.coordinates import SkyCoord
from collections import defaultdict

def setup_logging_to_file(filename):
	logging.basicConfig( filename='./'+filename,
						 filemode='w',
						 level=logging.INFO,
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

def build_filtered_table(table,flux, filter, filter_indices):
	if filter == 'False':
		return Table([table.ra.deg,table.dec.deg,flux], names=('RA','DEC','total_flux'))
	else: 
		df = Table([table.ra.deg,table.dec.deg,flux], names=('RA','DEC','total_flux'))
		RA = []
		DEC = []
		flux = []
		for i in range(len(filter_indices[0])):
			RA.append(np.average(df[filter_indices[0][i]]['RA']))
			DEC.append(np.average(df[filter_indices[0][i]]['DEC']))
			flux.append(np.average(df[filter_indices[0][i]]['total_flux']))
		return Table([RA,DEC,flux], names=('RA','DEC','total_flux'))

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
	pixcrd = np.array([[-10, -10], [24, 38], [45, 98]], np.float64)

	# Convert pixel coordinates to world coordinates
	world = w.wcs_pix2world(pixcrd, 1)

	# Convert the same coordinates back to pixel coordinates.
	pixcrd2 = w.wcs_world2pix(world, 1)

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

def write_correlation_params(table,prefix,correlator):
	'''
	Function writes the correlation parameters in a $SOURCE vex format or v2d format so that
	the correlator can read the phase centres in easily.
	e.g. def COSMOS-DEEP;
     source_name = COSMOS-DEEP;
	*    this source had calibrator code:  
	*    Center of COSMOS-deep
	     ra = 10h00m25.0000000s; dec =  02d33'00.000000"; ref_coord_frame = J2000;
	*    ra = 09h57m49.6816121s; dec =  02d47'25.904175"; ref_coord_frame = B1950;
	*    ra = 10h01m24.8107239s; dec =  02d27'22.307892"; ref_coord_frame = Date;
	enddef;
	Also returns the list of source names for inclusion in the $SCHED portion of the vex file.
	e.g., 
	scan No0066;
	*     Note a COMMENT was inserted during scheduling: 
	*       Loop 3 - part 1
     start=2022y265d00h08m26s; mode=EFF_BAND_32; source=R1_D;
     source=EK051E01;

	or for v2d:

	SOURCE A4038
	{
	doPointingCentre = True
	addPhaseCentre = name@A4038_A/RA@05:28:44.9836/Dec@-65:26:52.447
	addPhaseCentre = name@A4038_B/RA@05:28:44.6466/Dec@-65:26:44.711
	}
	'''
	if correlator == 'sfxc':
		correlation_string = ['$SOURCE;']
		source_string = ['$SCHED;']
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
			source_string.append('source=%s%s;'%(prefix[0:(8-sig_fig)],'{0:0{width}}'.format(i, width=sig_fig)))
		
		with open('%s_correlation_params.vex' % prefix, 'w') as f:
			for item in correlation_string:
				f.write("%s\n" % item)
			for item in source_string:
				f.write("%s\n" % item)
		f.close()
	if correlator == 'difx':
		correlation_string = ['SOURCE %s'%prefix,'{','doPointingCentre = True']
		for i in range(len(table['RA'])):
			sig_fig = len(str(len(table['RA'])))
			c = SkyCoord(table['RA'][i],table['DEC'][i],unit=('deg','deg'))
			sky_string = c.to_string('hmsdms').split(' ')
			RA_string = sky_string[0].replace('h',":").replace('m',":").replace('s','')
			Dec_string = sky_string[1].replace('d',":").replace('m',":").replace('s','')
			correlation_string.append('addPhaseCentre = name@%s%s/RA@%s/Dec@%s'%(prefix[0:(8-sig_fig)],'{0:0{width}}'.format(i, width=sig_fig),RA_string,Dec_string))
		correlation_string.append('}')
		with open('%s_correlation_params.v2d' % prefix, 'w') as f:
			for item in correlation_string:
				f.write("%s\n" % item)
		f.close()
	return

def primary_beam_power(offset, diameters, frequency, weights=None, pb_coeff=1.0):
	'''
	Analytic array primary-beam power response at an angular offset from the
	pointing centre, assuming a Gaussian voltage (element) beam per antenna.

	The voltage beam of antenna p is modelled as
	    E_p(r) = exp(-r^2 / 2 sigma_p^2),  sigma_p = FWHM_p / (2 sqrt(2 ln2)),
	    FWHM_p = pb_coeff * lambda / D_p.
	The baseline power beam is the (real) product
	    P_pq(r) = E_p E_q / sqrt(W_p W_q),
	which is itself a Gaussian with 1/sigma_pq^2 = 1/sigma_p^2 + 1/sigma_q^2.
	The array beam is the sensitivity-weighted sum over baselines
	    P_T(r) = sum_{j>i} w_ij exp(-r^2 / 2 sigma_ij^2),  w_ij = 1/sqrt(W_i W_j),
	and this returns the response normalised to the pointing centre, P_T(r)/P_T(0)
	(1.0 at centre, falling towards 0). Invert (P_T(0)/P_T(r)) for the correction.

	Parameters
	----------
	offset : float, array-like or astropy Quantity
	    Angular distance from the pointing centre. Plain numbers are degrees.
	diameters : float or array-like
	    Dish diameter(s) in metres. A scalar is treated as a homogeneous array
	    (the normalised shape is then independent of N_ant). A list/array gives a
	    heterogeneous array and the full per-baseline sum is evaluated.
	frequency : float or astropy Quantity
	    Observing frequency. Plain numbers are Hz.
	weights : array-like, optional
	    Per-antenna weights W_p (e.g. SEFD or noise variance). Only the relative
	    values affect the shape. Defaults to equal weighting. Ignored for the
	    homogeneous (scalar diameter) case.
	pb_coeff : float, optional
	    FWHM coefficient in FWHM = pb_coeff * lambda / D. 1.0 (default) for a
	    roughly uniform aperture; ~1.22 for a uniformly illuminated dish.

	Returns
	-------
	float or numpy.ndarray
	    Normalised primary-beam power (same shape as ``offset``), in [0, 1].
	'''
	# Coerce inputs to plain numbers in SI / radians
	if isinstance(offset, u.Quantity):
		r = offset.to(u.rad).value
	else:
		r = np.deg2rad(np.asarray(offset, dtype=float))
	if isinstance(frequency, u.Quantity):
		nu = frequency.to(u.Hz).value
	else:
		nu = float(frequency)

	c = 2.99792458e8           # speed of light, m/s
	lam = c / nu               # wavelength, m
	fwhm_to_sigma = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))

	D = np.atleast_1d(np.asarray(diameters, dtype=float))
	# Per-antenna Gaussian sigma (radians)
	sigma = pb_coeff * (lam / D) * fwhm_to_sigma

	# Homogeneous shortcut: all baselines identical, shape independent of N_ant
	if D.size == 1:
		sigma_bl2 = 0.5 * sigma[0] ** 2          # 1/sig_bl^2 = 2/sig^2
		return np.exp(-r ** 2 / (2.0 * sigma_bl2))

	# Heterogeneous array: explicit sum over unique baselines j > i
	n = D.size
	if weights is None:
		W = np.ones(n)
	else:
		W = np.atleast_1d(np.asarray(weights, dtype=float))
		if W.size != n:
			raise ValueError('weights must match the number of diameters')

	i_idx, j_idx = np.triu_indices(n, k=1)
	inv_sigma2 = 1.0 / sigma ** 2
	inv_sigma_bl2 = inv_sigma2[i_idx] + inv_sigma2[j_idx]   # 1/sigma_ij^2
	w_bl = 1.0 / np.sqrt(W[i_idx] * W[j_idx])               # baseline weight

	# Broadcast baselines against (possibly array-valued) offset
	r_arr = np.asarray(r, dtype=float)
	exponent = -0.5 * np.multiply.outer(r_arr ** 2, inv_sigma_bl2)
	P = np.sum(w_bl * np.exp(exponent), axis=-1) / np.sum(w_bl)
	return P if r_arr.ndim else float(P)

# ------------------------------------------------------------------------------
# Station look-up table + vex helpers for the array primary beam
# ------------------------------------------------------------------------------

# Single editable table of dish diameters / SEFDs, keyed by 2-letter station
# code. Regenerate from the upstream sources with build_stations.py.
STATIONS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
							 'data', 'stations.json')

def load_stations(path=STATIONS_FILE):
	'''Load the station look-up table (diameters + per-band SEFDs).

	Returns a dict keyed by upper-case 2-letter station code, e.g.
	    {'EF': {'diameter_m': 76.0, 'pb_model': 'G', 'sefd_jy': {'18': 19.0, ...}}}
	'''
	with open(path) as f:
		return json.load(f)

def _freq_to_band_cm(frequency, available):
	'''Return the entry from ``available`` (band labels in cm) whose wavelength
	is closest to ``frequency`` (Hz). ``available`` is an iterable of strings.'''
	lam_cm = 2.99792458e10 / float(frequency)   # c in cm/s over Hz
	bands = [b for b in available]
	return min(bands, key=lambda b: abs(float(b) - lam_cm))

def parse_vex_array(vexfile):
	'''Use vex.py to read which stations observe and the observing frequency.

	Returns ``(site_ids, frequency_hz)`` where ``site_ids`` is the list of
	2-letter station codes (in order of first appearance in the schedule).
	'''
	import io, contextlib
	sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
	from vex import Vex
	with contextlib.redirect_stdout(io.StringIO()):   # hush vex.py chatter
		v = Vex(vexfile)

	# Map site_name -> site_ID from the $SITE sector
	name2id = {}
	name = None
	for line in v.get_sector('SITE'):
		nm = v.get_variable('site_name', line)
		if nm:
			name = nm
		sid = v.get_variable('site_ID', line)
		if sid:
			name2id[name] = sid

	# Unique participating stations, in scan order
	seen = []
	for scan in v.sched:
		for s in scan['scan'].values():
			if s['site'] not in seen:
				seen.append(s['site'])

	site_ids = [name2id.get(nm, nm) for nm in seen]
	return site_ids, v.freq

def primary_beam_power_from_vex(offset, vexfile, frequency=None,
								weight_by_sefd=True, stations=None):
	'''Array primary-beam power at an angular ``offset`` for the array defined
	in a vex schedule.

	Pulls the participating stations (and observing frequency) from the vex file
	via vex.py, looks up each dish diameter and SEFD in the consolidated station
	table (data/stations.json), and evaluates :func:`primary_beam_power` with
	SEFD-weighted baselines.

	Parameters
	----------
	offset : float, array-like or astropy Quantity
	    Angular distance from the pointing centre (plain numbers are degrees).
	vexfile : str
	    Path to the .vex schedule.
	frequency : float or astropy Quantity, optional
	    Override the observing frequency (plain numbers are Hz). Defaults to the
	    value parsed from the vex file.
	weight_by_sefd : bool
	    If True (default) weight baselines by 1/sqrt(SEFD_i SEFD_j); falls back to
	    uniform weighting if any participating station lacks an SEFD.
	stations : dict, optional
	    Pre-loaded station table; defaults to :func:`load_stations`.

	Returns
	-------
	float or numpy.ndarray
	    Normalised primary-beam power (same shape as ``offset``), in [0, 1].
	'''
	if stations is None:
		stations = load_stations()

	site_ids, vex_freq = parse_vex_array(vexfile)
	if frequency is None:
		frequency = vex_freq
	elif isinstance(frequency, u.Quantity):
		frequency = frequency.to(u.Hz).value

	diameters, weights, missing = [], [], []
	for sid in site_ids:
		entry = stations.get(sid.upper())
		if entry is None or entry.get('diameter_m') is None:
			missing.append(sid)
			continue
		diameters.append(entry['diameter_m'])
		sefd = entry.get('sefd_jy') or {}
		if sefd:
			band = _freq_to_band_cm(frequency, sefd.keys())
			weights.append(sefd[band])
		else:
			weights.append(None)

	if missing:
		logging.warning('No dish diameter for station(s): %s - excluded from the '
						'primary beam' % ', '.join(missing))
	if len(diameters) < 2:
		raise ValueError('Need at least 2 stations with known diameters; '
						 'found %d' % len(diameters))

	if weight_by_sefd and all(w is not None for w in weights):
		W = weights
	else:
		if weight_by_sefd:
			logging.warning('Missing SEFD for some stations - using uniform weights')
		W = None

	return primary_beam_power(offset, diameters, frequency, weights=W)

# ------------------------------------------------------------------------------
# Thermal (radiometer-equation) sensitivity
# ------------------------------------------------------------------------------

def array_image_rms(sefds, bandwidth, t_int, eta=0.7, n_pol=2):
	'''Naturally-weighted image thermal noise of an interferometer (radiometer
	equation), combining all baselines:

	    dS_im = (1 / eta) * [ n_pol * sum_{j>i} 2 dnu t / (SEFD_i SEFD_j) ]^-1/2

	Parameters
	----------
	sefds : array-like
	    Per-antenna SEFDs (Jy).
	bandwidth : float or astropy Quantity
	    Total (synthesised) bandwidth; plain numbers are Hz.
	t_int : float or astropy Quantity
	    On-source integration time; plain numbers are seconds. Assumes every
	    baseline is present for this time (an upper bound for real schedules).
	eta : float
	    System/digitisation efficiency (~0.7 for 2-bit VLBI). Default 0.7.
	n_pol : int
	    Number of polarisations combined for Stokes I (2 gains sqrt(2)). Default 2.

	Returns
	-------
	float
	    1-sigma image rms in Jy/beam.
	'''
	if isinstance(bandwidth, u.Quantity):
		bandwidth = bandwidth.to(u.Hz).value
	if isinstance(t_int, u.Quantity):
		t_int = t_int.to(u.s).value

	S = np.atleast_1d(np.asarray(sefds, dtype=float))
	if S.size < 2:
		raise ValueError('Need at least 2 SEFDs')
	i_idx, j_idx = np.triu_indices(S.size, k=1)
	inv_var = np.sum(2.0 * bandwidth * t_int / (S[i_idx] * S[j_idx]))
	return 1.0 / (eta * np.sqrt(n_pol * inv_var))

def _scan_onsource_seconds(scan):
	'''On-source duration (s) of a single vex scan.'''
	starts = [s['scan_sec_start'] for s in scan['scan'].values()]
	ends = [s['scan_sec'] for s in scan['scan'].values()]
	if starts and ends:
		return max(ends) - min(starts)
	return 0.0


def _vex_onsource_seconds(sched, source=None, mk5clip=False,
						  mk5clip_seconds=1800.0):
	'''Total on-source time (s) summed over scans, optionally for one source.

	With ``mk5clip`` True, a leading block of consecutive same-source scans at
	the very start of the schedule is dropped if its cumulative on-source time
	exceeds ``mk5clip_seconds`` (default 30 min). This removes the strong
	fringe-finder / Mark5-clipped calibrator block that opens many schedules so
	it does not contribute to a wide-field rms estimate.
	'''
	start = 0
	if mk5clip and sched:
		lead_source = sched[0].get('source')
		block_end = 0
		while block_end < len(sched) and sched[block_end].get('source') == lead_source:
			block_end += 1
		block_time = sum(_scan_onsource_seconds(s) for s in sched[:block_end])
		if block_time > mk5clip_seconds:
			start = block_end

	total = 0.0
	for scan in sched[start:]:
		if source is not None and scan.get('source') != source:
			continue
		total += _scan_onsource_seconds(scan)
	return total

def expected_rms_from_vex(vexfile, bandwidth=None, source=None, t_int=None,
						  frequency=None, offset=0.0, eta=0.7, n_pol=2,
						  stations=None, mk5clip=False):
	'''Expected image thermal noise for the array in a vex schedule.

	Pulls the participating stations and (optionally) on-source time from the
	vex file, looks up each SEFD in the consolidated station table, and applies
	the radiometer equation via :func:`array_image_rms`. With ``offset`` > 0 the
	result is divided by the primary-beam power at that offset to give the
	effective noise away from the pointing centre.

	Parameters
	----------
	vexfile : str
	    Path to the .vex schedule.
	bandwidth : float or astropy Quantity, optional
	    Total synthesised bandwidth (plain numbers are Hz). If None it is taken
	    from the vex file (sum of the per-channel widths over unique sky-frequency
	    channels, i.e. ``Vex.total_bw_hz``).
	source : str, optional
	    Restrict the on-source time to this vex source name.
	t_int : float or astropy Quantity, optional
	    On-source integration time; if None it is summed from the schedule.
	frequency : float or astropy Quantity, optional
	    Override the observing frequency used for SEFD band selection.
	offset : float, array-like or astropy Quantity
	    Angular offset from the pointing centre for the effective noise (plain
	    numbers are degrees). Default 0.0 = central rms.
	eta, n_pol :
	    Passed through to :func:`array_image_rms`.
	stations : dict, optional
	    Pre-loaded station table; defaults to :func:`load_stations`.
	mk5clip : bool, optional
	    If True (and ``t_int`` is summed from the schedule), drop the leading
	    block of consecutive same-source scans at the very start of the
	    schedule when its cumulative on-source time exceeds 30 min. This
	    removes the strong fringe-finder / Mark5-clipped calibrator block that
	    opens many schedules so it does not contribute to the wide-field rms.
	    Default False.

	Returns
	-------
	float or numpy.ndarray
	    1-sigma rms in Jy/beam (scalar for central; array if ``offset`` is array).
	'''
	if stations is None:
		stations = load_stations()

	import io, contextlib
	sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
	from vex import Vex
	with contextlib.redirect_stdout(io.StringIO()):
		v = Vex(vexfile)

	# site_name -> site_ID, and participating stations in scan order
	name2id, name = {}, None
	for line in v.get_sector('SITE'):
		nm = v.get_variable('site_name', line)
		if nm:
			name = nm
		sid = v.get_variable('site_ID', line)
		if sid:
			name2id[name] = sid
	seen = []
	for scan in v.sched:
		for s in scan['scan'].values():
			if s['site'] not in seen:
				seen.append(s['site'])
	site_ids = [name2id.get(nm, nm) for nm in seen]

	if frequency is None:
		frequency = v.freq
	elif isinstance(frequency, u.Quantity):
		frequency = frequency.to(u.Hz).value
	if t_int is None:
		t_int = _vex_onsource_seconds(v.sched, source, mk5clip=mk5clip)
	if bandwidth is None:
		bandwidth = getattr(v, 'total_bw_hz', None)
		if bandwidth is None:
			raise ValueError('No total bandwidth in vex file; pass bandwidth '
							 'explicitly')

	sefds, missing = [], []
	for sid in site_ids:
		entry = stations.get(sid.upper())
		sefd = (entry or {}).get('sefd_jy') or {}
		if not sefd:
			missing.append(sid)
			continue
		band = _freq_to_band_cm(frequency, sefd.keys())
		sefds.append(sefd[band])
	if missing:
		logging.warning('No SEFD for station(s): %s - excluded from the rms'
						% ', '.join(missing))
	if len(sefds) < 2:
		raise ValueError('Need at least 2 stations with known SEFDs; found %d'
						 % len(sefds))

	rms = array_image_rms(sefds, bandwidth, t_int, eta=eta, n_pol=n_pol)

	if np.any(np.asarray(offset) != 0.0):
		rms = rms / primary_beam_power_from_vex(offset, vexfile, frequency=frequency,
												stations=stations)
	return rms

def locate_sources(vexfile):
    """Identify fringe finders, phase calibrators and target sources from a VEX schedule.

    Returns a tuple ``(fringe_finders, phase_refs, Sources)``.
    """
    fringe_finders = []
    phase_refs = []
    Sources = []

    ##fringe finders will be next to each other
    for i in range(len(vexfile.sched)-3): ### minus three to remove final scan which is two phase ref scans
        if vexfile.sched[i]['source']==vexfile.sched[i+1]['source'] and vexfile.sched[i]['scan'][0]['scan_sec']<300:
            fringe_finders.append(vexfile.sched[i]['source'])

    fringe_finders = list(set(fringe_finders))

    neighbors = defaultdict(set)
    time = defaultdict(float)

    ### phase reference and sources are one after another
    for i in range(len(vexfile.sched) - 1):
        s1 = vexfile.sched[i]['source']
        s2 = vexfile.sched[i+1]['source']
        t1 = vexfile.sched[i]['scan'][0]['scan_sec']
        t2 = vexfile.sched[i+1]['scan'][0]['scan_sec']

        time[s1]=t1
        time[s2]=t2

        if s1 != s2 and s1 not in fringe_finders and s2 not in fringe_finders:
            neighbors[s1].add(s2)
            neighbors[s2].add(s1)

        pairs = {tuple(sorted([s, list(neigh)[0]])) for s, neigh in neighbors.items() if len(neigh) == 1}

        pairs_with_times = {}

        for a, b in pairs:
            pairs_with_times[(a, b)] = (time[a], time[b])

    items = list(pairs_with_times.items())

    ##source will be longer observed than phase reference sources
    for i in range(len(items)):
        sources = items[i][0]
        times = items[i][1]
        if times[0]>times[1]:
            phase_refs.append(sources[1])
            Sources.append(sources[0])
        else:
            phase_refs.append(sources[0])
            Sources.append(sources[1])

    return fringe_finders, phase_refs, Sources