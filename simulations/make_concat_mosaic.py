import numpy as np
import os
import casatools
from casatasks import *
casalog.showconsole(True)
ms = []
for i in os.listdir('mosaic_ms'):
	print('Shifting %s'%i)
	os.system('rm -r mosaic_ms/%s_shift.ms'%(i.split('.ms')[0]))
	phaseshift(vis="mosaic_ms/%s"%i,outputvis='mosaic_ms/%s_shift.ms'%(i.split('.ms')[0]),phasecenter='J2000 12h00m00s +60d00m00')
	ms.append('mosaic_ms/%s_shift.ms'%(i.split('ms')[0]))
print('Concatenating')
concat(vis=ms,outputvis='mosaic_concat.ms')
np.save('mosaic_ms.npy',np.array(ms))
