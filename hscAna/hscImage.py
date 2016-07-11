"""
A collection of functions for getting 
and interacting with HSC images.
"""

import numpy as np
import pipeTools as pT

def write_fits(tract, patch, band='I', butler=None, save=True, outdir='/home/jgreco/projects/hscAna/output/'):
    if butler is None:
        butler = pT.get_butler()
    calexp = pT.get_exp(tract, patch, band, butler)
    outfile = outdir+'deepCoadd_HSC-'+band+'_'+str(tract)+'_'+patch+'.fits'
    calexp.writeFits(outfile)

def radec_to_tractpatch(ra, dec, butler=None, patch_as_str=True):
    """
    Get the tract and patch associated with the given ra and dec.

    Parameters
    ----------
    ra : float
        Right ascension of desired tract and patch.
    dec : float
        Declination of the desired tract and patch.
    butler : Bulter object, optional
        If None, create a butler object within function. 
        Otherwise, must be a butler object.
    patch_as_str : bool, optional
        If True, return patch as a string (e.g., '0,1').
        Otherwise, return patch as tuple (defualt = True).

    Returns
    -------
    tract : int
        HSC tract
    patch : string or tuple
        HSC patch
    """
    import lsst.afw.coord as afwCoord
    import lsst.afw.geom as afwGeom
    if butler is None:
        butler = pT.get_butler()
    skymap = butler.get('deepCoadd_skyMap', immediate=True)
    coord = afwCoord.IcrsCoord(afwGeom.Angle(ra, afwGeom.degrees), afwGeom.Angle(dec, afwGeom.degrees))
    tractInfo = skymap.findTract(coord)
    patchInfo = tractInfo.findPatch(coord)
    tract = tractInfo.getId()
    patch = patchInfo.getIndex()
    if patch_as_str:
        patch = str(patch[0])+','+str(patch[1])
    print tract, patch

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Get and save image given its tract and patch')
    parser.add_argument('tract', type=int, help='tract of observation')
    parser.add_argument('patch', type=str, help='patch of observation')
    parser.add_argument('-b', '--band', help='observation band', default='I')
    parser.add_argument('-o', '--outdir', help='output directory', 
                        default='/home/jgreco/projects/hscAna/output/')
    args = parser.parse_args()
    write_fits(args.tract, args.patch, band=args.band, outdir=args.outdir)
