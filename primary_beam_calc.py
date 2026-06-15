#!/usr/bin/env python3
'''
Standalone primary-beam attenuation calculator.

Give it a set of stations, an observing frequency, an (optional) bandwidth and
an angular offset from the pointing centre, and it returns the array primary-beam
power attenuation. Use --hpbw to also print the estimated array half-power beam
width (FWHM).

Examples
--------
    # EVN at 5 GHz, 12 arcmin off axis
    ./primary_beam_calc.py -s EF,JB,WB,ON,MC -f 5GHz -o 12arcmin --hpbw

    # Several offsets at once, band-averaged over 256 MHz
    ./primary_beam_calc.py -s EF,JB,WB -f 1.4GHz -b 256MHz -o 0,5,10,15arcmin

Stations, diameters and SEFDs are read from data/stations.json (the same table
used by the scheduler). Run with --list to see the available station codes.
'''
import argparse
import sys

import astropy.units as u

from wf_vlbi_functions import (load_stations, array_hpbw,
								primary_beam_power_from_stations)


def _quantity(value, default_unit):
	'''Parse a string like "5GHz" / "12 arcmin" / "1400" into an astropy
	Quantity, falling back to ``default_unit`` for a bare number.'''
	value = value.strip()
	try:
		return float(value) * default_unit
	except ValueError:
		return u.Quantity(value).to(default_unit)


def _parse_offsets(text):
	'''Parse a comma-separated list of offsets into a list of Quantities (deg).

	A unit on the final element applies to bare numbers earlier in the list, so
	"0,5,10,15arcmin" is read as four arcmin values.'''
	parts = [p.strip() for p in text.split(',') if p.strip()]
	# Find a trailing unit to apply to unitless entries.
	default = u.deg
	for p in reversed(parts):
		try:
			float(p)
		except ValueError:
			default = u.Quantity(p).unit
			break
	return [_quantity(p, default).to(u.deg) for p in parts]


def main(argv=None):
	ap = argparse.ArgumentParser(
		description='Primary-beam attenuation for a VLBI array.',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	ap.add_argument('-s', '--stations',
					help='Comma-separated station codes, e.g. EF,JB,WB,ON')
	ap.add_argument('-f', '--frequency',
					help='Centre frequency, e.g. 5GHz or 1400 (MHz assumed for '
						 'bare numbers)')
	ap.add_argument('-o', '--offset', default='0',
					help='Offset(s) from pointing centre, comma-separated, e.g. '
						 '12arcmin or 0,5,10arcmin (deg assumed for bare numbers)')
	ap.add_argument('-b', '--bandwidth', default=None,
					help='Total bandwidth for band-averaging, e.g. 256MHz '
						 '(MHz assumed for bare numbers). Omit for monochromatic.')
	ap.add_argument('--no-sefd-weight', action='store_true',
					help='Use uniform baseline weights instead of 1/sqrt(SEFD).')
	ap.add_argument('--pb-coeff', type=float, default=1.0,
					help='FWHM coefficient in FWHM = c * lambda / D.')
	ap.add_argument('--hpbw', action='store_true',
					help='Also report the estimated array half-power beam width.')
	ap.add_argument('--list', action='store_true',
					help='List the available station codes and exit.')
	args = ap.parse_args(argv)

	stations = load_stations()

	if args.list:
		print('Available stations (code: diameter):')
		for code in sorted(stations):
			d = stations[code].get('diameter_m')
			print('  %-4s %s m' % (code, d if d is not None else '?'))
		return 0

	if not args.stations or not args.frequency:
		ap.error('--stations and --frequency are required (or use --list)')

	frequency = _quantity(args.frequency, u.MHz)
	bandwidth = _quantity(args.bandwidth, u.MHz) if args.bandwidth else None
	offsets = _parse_offsets(args.offset)

	power, hpbw = primary_beam_power_from_stations(
		[o.to(u.deg).value for o in offsets], args.stations, frequency,
		bandwidth=bandwidth, weight_by_sefd=not args.no_sefd_weight,
		stations=stations, pb_coeff=args.pb_coeff, return_hpbw=True)

	codes = [c.strip().upper() for c in args.stations.split(',') if c.strip()]
	print('Stations : %s' % ', '.join(codes))
	print('Frequency: %s' % frequency.to(u.GHz))
	if bandwidth is not None:
		print('Bandwidth: %s (band-averaged)' % bandwidth.to(u.MHz))
	print('Weights  : %s' % ('uniform' if args.no_sefd_weight else '1/sqrt(SEFD)'))
	print()

	power = [power] if not hasattr(power, '__len__') else list(power)
	print('  %-14s %-12s %s' % ('offset', 'attenuation', '[mag]'))
	for off, p in zip(offsets, power):
		mag = '   (1/%.2f)' % (1.0 / p) if p > 0 else ''
		print('  %-14s %-12.4f%s' % (off.to(u.arcmin), p, mag))

	if args.hpbw:
		print()
		print('Estimated HPBW: %s' % hpbw.to(u.arcmin))

	return 0


if __name__ == '__main__':
	sys.exit(main())
