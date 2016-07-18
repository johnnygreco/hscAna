#!/usr/bin/env python 

"""
Get HSC deepCoadd fits files in a tract and patch that
contains a UDG candidate. 
"""

from __future__ import print_function

import os, sys, shutil
import numpy as np
import hscAna as ha
from lsst.daf.persistence import Butler
butler = Butler(ha.dataDIR)
band = 'I'

coords = np.loadtxt('../input/udg_candies.txt', skiprows=1, usecols=(0,1))

copydir = sys.argv[1]
outdir = '/home/jgreco/projects/hscAna/output/deepCoadds'

cmd = 'rsync -avc '+outdir+' '+copydir

for (ra, dec) in coords:
    tract, patch = ha.radec_to_tractpatch(ra, dec, butler=butler)

    print('getting deepCoadds for:', 'HSC-'+band+':', tract, patch)
    ha.write_deepCoadd_fits(tract, patch, band, butler=butler)

    # rsync fits files to different machine due to limited disk space
    os.system(cmd)

    # *carefully* delete the output directory
    if (outdir.split('/')[-2]=='output') and (outdir.split('/')[-1]=='deepCoadds'):
        print('deleting', outdir)
        shutil.rmtree(outdir)
    else:
        print('Careful!\n The output directory is not what you think!')
