from utils import *
from cattools import *
from cuts import *

tract = 9347
patch = '5,7'
band = 'I'
butler = get_butler()
cat = butler.get('deepCoadd_meas', tract=tract, patch=patch, filter='HSC-'+band, immediate=True)
exp = butler.get('deepCoadd_calexp', tract=tract, patch=patch, filter='HSC-'+band, immediate=True)
wcs = exp.getWcs()
calib = exp.getCalib()
cut = basic_cut(cat)



