import numpy as np
from astropy.table import Table
from astroquery.utils.tap.core import TapPlus
import astropy.units as u
import sys
import json

with open(sys.argv[1]) as f:
	inputs = json.load(f)

ra = inputs['ra']
dec = inputs['dec']
FoV = inputs['rad_FoV']

casdatap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")

job = casdatap.launch_job_async("SELECT * FROM AS110.racs_mid_components_v01 where 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{},{},{}))".format(ra,dec,FoV))

r = job.get_results()

keep_cols = ['ra','dec','total_flux']
r=r[keep_cols]
print(r)
r.write('RACS_test.csv',overwrite=True)

