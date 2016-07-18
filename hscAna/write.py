
from __future__ import division, print_function

__all__ = ['make_default_outdir', 'write_deepCoadd_fits']

def make_default_outdir(tract, patch, band):
    """
    Make the default output directory for the deepCoadd
    fits files that are the output of write_deepCoadd_fits. 
    The default directory has the following structure:
    ../output/deepCoadds/HSC-band/tract/patch.
    """
    import os

    # output directory is in previous directory 
    outdir = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.dirname(outdir)

    dirs = ['output', 'deepCoadds', 'HSC-'+band.upper(),
            str(tract), patch[0]+'-'+patch[-1]]

    # build path, make directories if they don't exist
    for d in dirs:
        outdir = os.path.join(outdir, d)
        if not os.path.isdir(outdir):
            print('created', outdir)
            os.mkdir(outdir)

    return outdir

def write_deepCoadd_fits(tract, patch, band='I', outdir='default', butler=None):
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
    outdir : string, optional
        The output directory. If 'default', the directory 
        will be the output of make_default_outdir.
    butler : Butler object, optional
        If None, a Butler object will be created. 
    
    Notes
    -----
    The output files will be written to outdir/folder, with names
    1) img.fits (image)
    2) bad.fits (bad pixel mask)
    3) det.fits (detection pixel mask)
    4) sig.fits (sigma image)
    5) psf.fits (point spread function)
    """
    import os
    from astropy.io import fits
    from myPipe import MyPipe

    band = band.upper()
    pipe = MyPipe(tract, patch, band=band, butler=butler)

    if outdir=='default':
        outdir = make_default_outdir(tract, patch, band)

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

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Get and save deepCoadd image in given tract and patch')
    parser.add_argument('tract', type=int, help='tract of observation')
    parser.add_argument('patch', type=str, help='patch of observation')
    parser.add_argument('-b', '--band', help='observation band', default='I')
    parser.add_argument('-o', '--outdir', help='output directory', default='default')
    args = parser.parse_args()
    write_deepCoadd_fits(args.tract, args.patch, band=args.band, outdir=args.outdir)
