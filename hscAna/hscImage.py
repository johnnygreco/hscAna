"""
A collection of functions for getting 
and interacting with HSC images.
"""

from __future__ import division, print_function

__all__ = ['write_deepCoadd_fits', 'radec_to_tractpatch']

def write_deepCoadd_fits(tract, patch, band='I', folder=None, outdir='../output'):
    """
    Write deepCoadd fits images for the given tract, patch, and band.
    Will write individual files for the image, bad pixel mask, detected
    pixel mask, and sigma map. 

    Parameters
    ----------
    tract : int
        HSC tract.
    patch : sting
        HSC patch. 
    band : string, optional
        HSC filter (GRIZY).
    folder : string, optional
        Name of folder for output files. If None, one will be
        created called deepCoadd_tract_patch_HSC-band.
    outdir : string, optional
        The output directory. 
    
    Notes
    -----
    The output files will be written to outdir/folder, with names
    1) img.fits (image)
    2) bad.fits (bad pixel mask)
    3) det.fits (detection pixel mask)
    4) sig.fits (sigma image)
    """
    import os
    from astropy.io import fits
    from myPipe import MyPipe

    band = band.upper()
    pipe = MyPipe(tract, patch, band=band)

    if folder is None:
        folder = 'deepCoadd_{}_{}_HSC-{}'.format(tract, patch, band)
    outdir = os.path.join(outdir, folder)

    # make output directory if it doesn't exist
    if not os.path.isdir(outdir):
        print('created', outdir)
        os.mkdir(outdir)

    # get headers: 0=image, 1=mask, 2=variance
    hdulist= fits.open(pipe.get_fn())
    headers = [hdulist[i].header for i in range(1,4)]
    headers[0].set('PIXSCALE', pipe.get_pixscale())
    headers[0].set('ZP_PHOT', pipe.get_zptmag())

    # write images
    getters = [pipe.get_img, pipe.get_badmask, pipe.get_detmask, pipe.get_sigma]
    labels = ['img', 'bad', 'det', 'sig']
    headnums = [0, 1, 1, 2]
    for get, lab, num in zip(getters, labels, headnums):
        print('writing', lab+'.fits')
        fn = os.path.join(outdir, lab+'.fits')
        fits.writeto(fn, get(), headers[num], clobber=True)

    # write psf fits file
    fn = os.path.join(outdir, 'psf.fits')
    print('writing', 'psf.fits')
    pipe.calexp.getPsf().computeImage().writeFits(fn)


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
        import lsst.daf.persistence
        from myPipe import dataDIR
        butler = lsst.daf.persistence.Butler(dataDIR)
    skymap = butler.get('deepCoadd_skyMap', immediate=True)
    coord = afwCoord.IcrsCoord(afwGeom.Angle(ra, afwGeom.degrees), afwGeom.Angle(dec, afwGeom.degrees))
    tractInfo = skymap.findTract(coord)
    patchInfo = tractInfo.findPatch(coord)
    tract = tractInfo.getId()
    patch = patchInfo.getIndex()
    if patch_as_str:
        patch = str(patch[0])+','+str(patch[1])
    print(tract, patch)
    return tract, patch

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Get and save deepCoadd image in given tract and patch')
    parser.add_argument('tract', type=int, help='tract of observation')
    parser.add_argument('patch', type=str, help='patch of observation')
    parser.add_argument('-b', '--band', help='observation band', default='I')
    parser.add_argument('-o', '--outdir', help='output directory', 
                        default='/home/jgreco/projects/hscAna/output/')
    args = parser.parse_args()
    write_deepCoadd_fits(args.tract, args.patch, band=args.band, outdir=args.outdir)
