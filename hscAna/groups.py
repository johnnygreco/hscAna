"""
Some tools to get HSC data associated 
with galaxy groups. 
"""

from __future__ import print_function

__all__ = ['get_group_fits']

import numpy as np

def get_group_fits(ra, dec, z, group_id, box_width=3.0, band='I', outdir=None, butler=None):
    """
    Get fits files within width/2 of the given coords.  

    Parameters
    ----------
    ra, dec, z : float
        The luminosity-weighted right ascension, 
        declination and redshift of a galaxy group.
    group_id : int
        The galaxy group id. 
    box_width : float, optional 
        The width of the data region in Mpc.
    band : string, optional
        The photometric band (GRIZY). 
    outdir : string, optional
        If not None, the output directory, else 
        it will be ../output/
    butler : Butler object
        If None, a butler will be created.

    Notes
    -----
    If it does not exist, a directory called group_id will 
    be created in outdir, and he fits files will be saved
    there. 
    """
    import os, shutil
    import utils
    from myPipe import dataDIR
    from params.copydir import copydir
    from write import write_deepCoadd_fits
    from toolbox.cosmo import Cosmology

    if butler is None:
        import lsst.daf.persistence
        butler = lsst.daf.persistence.Butler(dataDIR)
    if outdir is None: 
        outdir = os.path.dirname(os.path.abspath(__file__))
        outdir = os.path.join(os.path.dirname(outdir), 'output')
    outdir = os.path.join(outdir, 'group_'+str(group_id))
    if not os.path.isdir(outdir):
        print('creating', outdir)
        os.mkdir(outdir)

    cmd = 'rsync -avc '+outdir+' '+copydir

    D_A = Cosmology().D_A(z) # angular diameter distance
    theta = (box_width/D_A)*180.0/np.pi
    print('will extract a sky box with sides of ', theta, 'degrees')
    regions = utils.get_hsc_regions(utils.skybox(ra, dec, theta), butler=butler)
    print('***** found', len(regions), 'frames in region *****')

    for tract, patch in regions:
        print('getting deepCoadds for:', 'HSC-'+band+':', tract, patch)
        prefix = 'HSC-'+band.upper()+'_'+str(tract)+'_'+patch[0]+'-'+patch[-1]

        write_deepCoadd_fits(tract, patch, band, butler=butler, outdir=outdir, prefix=prefix)

        # rsync fits files to different machine due to limited disk space
        os.system(cmd)

        # delete files
        files = os.listdir(outdir)
        for f in files:
            print('deleting', f)
            os.remove(os.path.join(outdir, f))
    print('deleting', outdir)
    shutil.rmtree(outdir)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Get and deepCoadd fits files near a galaxy group')
    parser.add_argument('ra', type=float, help='luminosity-weighted ra of group')
    parser.add_argument('dec', type=float, help='luminosity-weighted dec of group')
    parser.add_argument('z', type=float, help='luminosity-weighted redshift of group')
    parser.add_argument('group_id', type=int, help='group id')
    parser.add_argument('-w', '--box_width', help='width of the data region in Mpc', default=3.0)
    parser.add_argument('-b', '--band', help='observation band', default='I')
    parser.add_argument('-o', '--outdir', help='output directory', default=None)
    args = parser.parse_args()
    get_group_fits(args.ra, args.dec, args.z, args.group_id, args.box_width, args.band, args.outdir)
