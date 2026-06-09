#!/usr/bin/env python3
"""Regenerate data/stations.json - the single editable table of dish diameters
and per-band SEFDs used by wf_vlbi_functions.primary_beam_power*.

Diameters / beam models come from the VLBI_pipeline repo and SEFDs from the
antab_editor repo. Hand-maintained additions (e.g. LBA dishes absent upstream)
live in LBA_DIAMETERS below. After running, edit data/stations.json directly for
any local tweaks.
"""
import json
import os
from urllib.request import urlopen

PB_URL   = 'https://raw.githubusercontent.com/jradcliffe5/VLBI_pipeline/master/data/primary_beams.json'
SEFD_URL = 'https://raw.githubusercontent.com/jradcliffe5/antab_editor/main/sefd_values.txt'

# Codes appearing in SEFD table that map onto a different diameter-table code
SEFD_ALIAS = {'JB1': 'JB', 'JB2': 'JB2'}

# Dishes missing from the EVN-centric diameter table (metres)
LBA_DIAMETERS = {'AT': 22.0, 'PA': 64.0, 'MP': 22.0, 'CD': 30.0, 'HO': 26.0}

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'stations.json')


def main():
    pb = json.loads(urlopen(PB_URL).read().decode())
    stations = {}
    for code, bands in pb.items():
        band = bands.get('L', next(iter(bands.values())))
        stations[code.upper()] = {
            'diameter_m': band['diameter'],
            'pb_model':   band.get('pb_model', 'G'),
            'sefd_jy':    {},
        }

    sefd_txt = urlopen(SEFD_URL).read().decode()
    rows = [r for r in sefd_txt.splitlines() if r.strip()]
    bands_cm = [c.strip() for c in rows[0].split('|')[1:]]
    for row in rows[1:]:
        cells = [c.strip() for c in row.split('|')]
        code = cells[0].upper()
        code = SEFD_ALIAS.get(code, code)
        sefd = {b: float(v) for b, v in zip(bands_cm, cells[1:]) if v}
        if not sefd:
            continue
        stations.setdefault(code, {'diameter_m': None, 'pb_model': 'G', 'sefd_jy': {}})
        stations[code]['sefd_jy'].update(sefd)

    for code, d in LBA_DIAMETERS.items():
        stations.setdefault(code, {'diameter_m': None, 'pb_model': 'G', 'sefd_jy': {}})
        stations[code]['diameter_m'] = d

    stations = dict(sorted(stations.items()))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w') as f:
        json.dump(stations, f, indent=2)

    nd = sum(1 for s in stations.values() if s['diameter_m'])
    ns = sum(1 for s in stations.values() if s['sefd_jy'])
    print('wrote %s: %d stations (%d with diameter, %d with SEFD)'
          % (OUT, len(stations), nd, ns))


if __name__ == '__main__':
    main()
