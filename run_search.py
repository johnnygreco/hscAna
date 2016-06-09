#!/usr/bin/env python 

import sys
import hscana
import numpy as np 
import argparse
from toolbox.astro import angsep
from astropy.table import Table
from hscana.utils import get_hsc_regions, skybox
group_info = Table.read('/home/jgreco/data/groups/group_info.csv')

def run(group_id, band='I', box_width=3.0, max_sep=2.0, butler=None):
    """
    Run script to search for UDG candidates.

    Parameters
    ----------
    group_id : int
        Galaxy group identification number.
    band : string, optional, default is 'I'
        HSC band (G, R, I, Z, or Y).
    box_width : float, optional, default is 3 degree
        Angular width of 'skybox' to extract around the groups
        luminosity weighted center. Must be in degrees.
    max_sep : float, optional, default is 2 arcsec
        Maximum separation for which two objects are considered
        to be the same object. Must be in arcseconds. 
    butler : Butler object, optional
        If None, will be created within function. If you are
        looping over many groups, you should create a butler
        once outside this function. 
    """
    if butler is None:
        butler = hscana.get_butler()
    idx = np.argwhere(group_info['group_id']==group_id)[0,0]
    ra_c, dec_c, D_A, D_L = group_info['ra', 'dec', 'D_A', 'D_L'][idx]
    theta = (box_width/D_A)*180.0/np.pi
    group_regions = get_hsc_regions(skybox(ra_c, dec_c, theta, theta))
    print 'We will extract region of angular size theta =', round(theta, 3), 'degree'

    # build list of udg candidates
    candy = []
    for tract, patch in group_regions:
        # some tracts and patches are missing
        print tract, patch
        try:
            mycat = hscana.MyCat(tract, patch, band, group_id, makecuts=True, butler=butler)
        except:
            print '!!!!! FAILED !!!!!'
            continue
        if mycat.count()>0:
            candy.append(mycat)
            
    # get coordinates of all candidates
    coords = []
    for i in range(len(candy)):
        for ra, dec in candy[i].coord():
            coords.append((ra, dec))
    coords = np.array(coords)
            
    # build mask for double entries
    mask = np.ones(coords.shape[0], dtype=bool)
    for i, (ra, dec) in enumerate(coords):
        if mask[i]==True:
            unique = angsep(ra, dec, coords[:,0], coords[:,1]) > max_sep
            unique[i] = True
            mask &= unique
    coords = coords[mask]
    print 'number of candidates =', coords.shape[0]

    # output in format for hscMap
    if coords.shape[0]>0:
        np.savetxt('output/group_'+str(group_id)+'_coords.csv', coords, delimiter=',', header='ra,dec')
    else:
        print 'group', group_id, 'has zero candidates'

if __name__=='__main__':
    # for usage, enter ./run_search.py
    parser = argparse.ArgumentParser(description='Search for UDG candidates.')
    parser.add_argument('-g', '--group_id', type=int, default=None, help='run search for single group with id=group_id')
    parser.add_argument('-n', '--Ngal', type=int, default=None, help='run search on all groups with <= Ngal galaxies')
    parser.add_argument('-z', '--z', type=float, default=None, help='run search on all groups with redshift < z')
    args = parser.parse_args()
    if args.group_id is not None:
        print 'running search for group', args.group_id
        run(group_id=args.group_id)
    elif args.Ngal is not None:
        butler = hscana.get_butler()
        cut = group_info['Ngal']<=args.Ngal
        print 'running search for all groups with Ngal <=', args.Ngal
        if args.z is not None:
            cut &= group_info['z']<args.z
            print 'and z <', args.z
        for ID in group_info[cut]['group_id']:
            print '***** searching in group '+str(ID)+' *****'
            run(group_id=ID, butler=butler)
    elif args.z is not None:
        butler = hscana.get_butler()
        cut = group_info['z']<args.z
        print 'running search for all groups with z <', args.z
        for ID in group_info[cut]['group_id']:
            print '***** searching in group '+str(ID)+' *****'
            run(group_id=ID, butler=butler)
    else:
        parser.print_help()
