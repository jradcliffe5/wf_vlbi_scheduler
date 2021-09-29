import numpy as np
import os, ast
import casatools
from casatasks import *
from simulator_functions import headless
casalog.showconsole(True)

inputs = headless('simulator_inputs.txt')

## get images
pbs = []
ims = []
for i in os.listdir('%s'%inputs['output_path']):
	if i.endswith('.ms.pb'):
		pbs.append("%s/%s"%(inputs['output_path'],i))
	elif i.endswith('.ms_pb.image'):
		ims.append("%s/%s"%(inputs['output_path'],i))

lm = casatools.linearmosaic()
os.system('rm -r test.linmos test.weightlinmos')
print('Defining output image')
imsize = np.array(ast.literal_eval(inputs['mosaic_area']))*1.2*3.6e6
max_mosaic_size = 7200
cell = '%.5farcsec'%(imsize/max_mosaic_size)[0]
pc = ast.literal_eval(inputs['field_centre'])
lm.defineoutputimage(nx=max_mosaic_size, cellx=cell, imagecenter='%s %s'%(pc[0],pc[1]), outputimage='%s/mosaic.linmos'%inputs['output_path'], outputweight='%s/mosaic.weightlinmos'%inputs['output_path'])
print('Making linear mosaic')
for i in range(len(ims)):
	print(i)
	lm.makemosaic(images=ims[i], weightimages=pbs[i],weighttype=2,imageweighttype=0)
exportfits(imagename='%s/mosaic.linmos'%inputs['output_path'],fitsimage='%s/mosaic.linmos.fits'%inputs['output_path'],overwrite=True)
exportfits(imagename='%s/mosaic.weightlinmos'%inputs['output_path'],fitsimage='%s/mosaic.weightlinmos.fits'%inputs['output_path'],overwrite=True)