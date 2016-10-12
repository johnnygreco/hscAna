
from __future__ import division, print_function

__all__ = ['skybox', 'get_hsc_regions', 'radec_to_tractpatch']

import numpy as np

def skybox(ra_c, dec_c, width, height=None):
    """
    Calculate the four corners of a box centered at (ra_c, dec_c).
    All input parameters must be in degrees.

    Parameters
    ----------
    ra_c : float
        The central ra of the box in degrees.
    dec_c : float
        The central dec of the box in degrees.
    width : float
        The angular width of the box in degrees.
    height : float, optional
        The angular height of the box in degrees.
        If None, will set height=width. 

    Returns
    -------
    box_coords : list of tuples
        The four coordinates for the box corners.

    Note
    ----
    This calculation is only an approximation, as it 
    assumes the angular separation is small. 
    """
    if height is None:
        height = width
    dec_lo = dec_c - height/2.0
    dec_hi = dec_c + height/2.0
    ra_min_lo = ra_c - width/2.0/np.cos(dec_lo*np.pi/180.)
    ra_min_hi = ra_c - width/2.0/np.cos(dec_hi*np.pi/180.)
    ra_max_lo = ra_c + width/2.0/np.cos(dec_lo*np.pi/180.)
    ra_max_hi = ra_c + width/2.0/np.cos(dec_hi*np.pi/180.)
    box_coords = [(ra_min_lo, dec_lo),
                  (ra_min_hi, dec_hi),
                  (ra_max_hi, dec_hi),
                  (ra_max_lo, dec_lo)]
    return box_coords

def get_hsc_regions(box_coords, butler=None):
    """
    Get hsc regions within a polygonal region (box) of the sky. Here, 
    hsc regions means the tracts and patches within the 'skybox'. 

    Parameters
    ----------
    box_coords : list of tuples
        The four coordinates for the box corners. This is the 
        output of the skybox function. If only one coordinate
        is given, will return a single tract and patch. 
    butler : Butler object, optional
        If None, then a will be created in this function.
        Default is None.

    Returns
    -------
    regions : structured ndarray
        The tracts and patches within the skybox. The columns of 
        the array are 'tract' and 'patch'.

    Note
    ----
    This may give incorrect answers on regions that are larger than a tract, 
    which is ~1.5 degree = 90 arcminute.
    """
    import lsst.afw.coord as afwCoord
    import lsst.afw.geom as afwGeom
    if butler is None:
        import lsst.daf.persistence
        from myPipe import dataDIR
        butler = lsst.daf.persistence.Butler(dataDIR)
    if len(box_coords)==4:
        from toolbox.astro import angsep
        (ra1, dec1), (ra2, dec2) = box_coords[0], box_coords[2]
        if angsep(ra1, dec1, ra2, dec2, sepunits='arcmin') > 90.0:
            print('\n********* WARNING *********')
            print('Region larger than a tract')
            print('***************************\n')
    skymap = butler.get('deepCoadd_skyMap', immediate=True)
    coordList = [afwCoord.IcrsCoord(afwGeom.Angle(ra, afwGeom.degrees),\
                 afwGeom.Angle(dec, afwGeom.degrees)) for ra, dec in box_coords]
    tractPatchList = skymap.findClosestTractPatchList(coordList)
    regions = []
    for tractInfo, patchInfoList in tractPatchList:
        for patchInfo in patchInfoList:
            patchIndex = patchInfo.getIndex()
            regions.append((tractInfo.getId(), str(patchIndex[0])+','+str(patchIndex[1])))
    return np.array(regions, dtype=[('tract', int), ('patch', 'S4')])

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
    return tract, patch
