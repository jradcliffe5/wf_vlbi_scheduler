import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from simulator_functions import headless
from astropy.wcs import WCS
from matplotlib.gridspec import GridSpec


inputs = headless('simulator_inputs.txt')


hdu = fits.open('single_pointing-beam-I.fits')
cdelt=np.abs(hdu[0].header['CDELT1'])
data = hdu[0].data.squeeze()
w = WCS(hdu[0].header,naxis=2)
y, x = np.mgrid[:data.shape[0], :data.shape[1]]


# In[17]:


# Fit the data using a Gaussian
g_init = models.Gaussian2D(amplitude=1., x_mean=data.shape[0]//2,y_mean=data.shape[1]//2, x_stddev=1.,y_stddev=1.)
g_init.amplitude.fixed = True
fit_g = fitting.LevMarLSQFitter()
g_air = models.AiryDisk2D(amplitude=1., x_0=data.shape[0]//2,y_0=data.shape[1]//2, radius=1)
g_air.amplitude.fixed = True
g = fit_g(g_init, x,y,data,weights=(data/data.max())**2)
g_air = fit_g(g_air, x,y,data,weights=(data/data.max())**2)


# In[18]:


fig = plt.figure(figsize=(27, 10))
gs = GridSpec(nrows=1,ncols=4,wspace=0.05)
ax = fig.add_subplot(gs[0],projection=w)
ax.imshow(data, origin='lower', interpolation='nearest', vmin=0, vmax=1,rasterized=True)
ax.set_title("Data")
ax.set_ylabel('Declination (J2000)')
ax.set_xlabel('')
ax = fig.add_subplot(gs[1],projection=w)
ax.imshow(g(x, y), origin='lower', interpolation='nearest', vmin=0, vmax=1,rasterized=True)
ax.set_yticklabels([])
ax.set_ylabel('')
ax.set_xlabel('Right Ascension (J2000)')
ax.set_title("Model (Gaus)")
ax = fig.add_subplot(gs[2],projection=w)
ax.imshow(g_air(x, y), origin='lower', interpolation='nearest', vmin=0, vmax=1,rasterized=True)
ax.set_title("Model (Airy)")
ax.set_yticklabels([])
ax.set_ylabel('')
ax.set_xlabel('')
ax = fig.add_subplot(gs[3],projection=w)
ax.imshow(data - g(x, y), origin='lower', interpolation='nearest',rasterized=True)
ax.set_title("Residual (Gaus)")
ax.set_yticklabels([])
ax.set_ylabel('')
ax.set_xlabel('')
fig.savefig('%s/Primary_beam_fit_2D.pdf'%(inputs['output_path']))
plt.clf()

# In[19]:


fig = plt.figure()
ax = fig.add_subplot(111)
_x = np.linspace(-1*cdelt*data.shape[0]/2,cdelt*data.shape[1]/2,data.shape[0])
ax.plot(_x, data[data.shape[0]//2,:],c='k',label='PB')
ax.plot(_x,g(x, y)[data.shape[0]//2,:],c='b',ls='--',label='Gaus')
ax.plot(_x,g_air(x, y)[data.shape[0]//2,:],c='r',ls=':',label='Airy')
ax.legend()
fig.savefig('%s/Primary_beam_fit_1D.pdf'%(inputs['output_path']))
plt.clf()


# In[25]:


fwhm = np.mean([g.x_stddev.value*2*np.sqrt(2*np.log(2))*cdelt,g.y_stddev.value*2*np.sqrt(2*np.log(2))*cdelt])
print(fwhm)
np.save("%s/PB_fit.npy"%(inputs['output_path']),fwhm)

