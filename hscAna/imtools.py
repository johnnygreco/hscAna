from __future__ import print_function

from astropy.io import fits

def sig_to_wts(sigfile, wfile='wts.fits'):
    """
    Convert sigma image to weights image, where 
    weight = 1/sigma**2. 

    Parameters
    ----------
    sigfile : string
        The input sigma image file.
    wfile : string, optional
        The output weights image file.
    """
    sigfits = fits.open(sigfile)[0]
    weights = 1.0/sigfits.data**2
    print('writing', wfile)
    fits.writeto(wfile, weights, sigfits.header, clobber=True)


def wts_with_badpix(wfile, badfile, wnewfile='wts_bad.fits', flagval=-100.0):
    """
    Flag bad pixels in the weight image for sextractor.

    Parameters
    ----------
    wfile : string
        Input weight image file (weight = 1/sigma**2).
    badfile : string
        Input bad pixel map file (0 = good pixels).
    wnewfile : string, optional
        Output weights + badpix file.
    flagval : float, optional
        The weight to be assigned to bad pixels.
        (sextractor default threshhold = 0)
    """
    badpix = fits.getdata(badfile)
    wfits = fits.open(wfile)[0]
    wfits.data[badpix!=0] = flagval
    print('writing', wnewfile)
    fits.writeto(wnewfile, wfits.data, wfits.header, clobber=True)
