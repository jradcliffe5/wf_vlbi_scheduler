import numpy as np
import os
import casatools
from casatasks import *
casalog.showconsole(True)
import sys

ms = sys.argv[1]

importfits(fitsimage='%s_image-beam-I.fits'%ms,imagename="%s.pb"%ms)
importfits(fitsimage='%s_image-image-pb.fits'%ms,imagename='%s_pb.image'%ms)