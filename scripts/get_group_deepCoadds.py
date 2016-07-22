#!/usr/bin/env python 

"""
Get HSC deepCoadd fits files for Yang galaxy groups.
"""

from __future__ import print_function

import os, sys, shutil
import numpy as np
from astropy.table import Table
import hscAna

import lsst.daf.persistence
butler = lsst.daf.persistence.Butler(hscAna.dataDIR)

group_info = Table.read('/home/jgreco/data/groups/group_info.csv')
print(len(group_info), 'galaxy groups in catalog')

##############################################################
# Make cuts on group catalog
##############################################################
zmax = 0.08
Ngal_max = 10
cut  = group_info['z'] <= zmax
cut &= group_info['Ngal'] <= Ngal_max
group_info = group_info[cut]
print(len(group_info), 'galaxy groups in catalog')

##############################################################
# Get tracts and patches for each group and associated files
##############################################################

band = 'I'
box_width = 3.0 # Mpc

copydir = sys.argv[1]
deepCoadds_dir = '/home/jgreco/projects/hscAna/output/deepCoadds'
cmd = 'rsync -avc '+deepCoadds_dir+' '+copydir

for group_id, ra, dec, D_A in group_info['group_id', 'ra','dec', 'D_A']:
    theta = (box_width/D_A)*180.0/np.pi
    group_regions = hscAna.get_hsc_regions(hscAna.skybox(ra, dec, theta), butler=butler)

    for tract, patch in group_regions:
        print('getting deepCoadds for:', 'HSC-'+band+':', tract, patch)
        outdir = deepCoadds_dir

        # make output directories if they don't exist
        # path will be {outdir}/HSC-band/group_id/tract/patch
        dirs = ['', 'HSC-'+band, 'group_'+str(group_id), 
                str(tract), patch[0]+'-'+patch[-1]]
        for d in dirs:
            outdir = os.path.join(outdir, d)
            if not os.path.isdir(outdir):
                print('created', outdir)
                os.mkdir(outdir)

        hscAna.write_deepCoadd_fits(tract, patch, band, butler=butler, outdir=outdir)

        # rsync fits files to different machine due to limited disk space
        os.system(cmd)

        # *carefully* delete the output directory
        if (deepCoadds_dir.split('/')[-2]=='output') and (deepCoadds_dir.split('/')[-1]=='deepCoadds'):
            print('deleting', deepCoadds_dir)
            shutil.rmtree(deepCoadds_dir)
        else:
            print('***** Careful! ***** \n The deepCoadds_dir is not what you think!')
