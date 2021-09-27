import numpy as np
import os
import casatools
from casatasks import *
casalog.showconsole(True)

## get images
pbs = []
ims = []
for i in os.listdir('mosaic_ms'):
	if i.endswith('.ms.pb'):
		pbs.append("mosaic_ms/%s"%i)
	elif i.endswith('.ms_pb.image'):
		ims.append("mosaic_ms/%s"%i)

lm = casatools.linearmosaic()
os.system('rm -r test.linmos test.weightlinmos')
print('Defining output image')
lm.defineoutputimage(nx=7200, cellx='1arcsec', imagecenter='12h00m00s 60d00m00s', outputimage='test.linmos', outputweight='test.weightlinmos')
print('Making linear mosaic')
for i in range(len(ims)):
	print(i)
	lm.makemosaic(images=ims[i], weightimages=pbs[i],weighttype=2,imageweighttype=0)
exportfits(imagename='test.linmos',fitsimage='test.linmos.fits',overwrite=True)
exportfits(imagename='test.weightlinmos',fitsimage='test.weightlinmos.fits',overwrite=True)