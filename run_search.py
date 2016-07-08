#!/usr/bin/env python 

import numpy as np
from astropy.table import Table
import hscAna
from hscAna.search import group_search
group_info = Table.read('/home/jgreco/data/groups/group_info.csv')
print len(group_info), 'galaxy groups in catalog'

gemini_cuts = True

##########################################################
# Make cuts on group catalog
##########################################################
if gemini_cuts:
    zmax = 0.08
    Ngal_max = 15
    # only look at ra's accessible to gemini
    cut  = (group_info['ra'] > 12.5*15) | (group_info['ra'] < 5.0*15)
    cut &= group_info['z'] <= zmax
    cut &= group_info['Ngal'] <= Ngal_max
    group_info = group_info[cut]
    print len(group_info), 'galaxy groups after cuts'

##########################################################
# Perform UDG search within all remaining groups
##########################################################
butler = hscana.get_butler()
for ID in group_info['group_id']:
    print '***** searching in group '+str(ID)+' *****'
    idx = np.argwhere(group_info['group_id']==ID)[0,0]
    ra, dec, z, Ngal = group_info['ra', 'dec', 'z', 'Ngal'][idx]
    print 'ra dec =', ra, dec
    print 'z =', z
    print 'Ngal=', Ngal 
    group_search(group_id=ID, butler=butler)
