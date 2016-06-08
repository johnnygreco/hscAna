import numpy as np

def skybox(ra_c, dec_c, width, height):
    """
    Calculate the four corners of a box centered at (ra_c, dec_c).
    All input parameters must be in degrees.

    Parameters
    ----------
    ra_c: float
        The central ra of the box in degrees.
    dec_c: float
        The central dec of the box in degrees.
    width: float
        The angular width of the box in degrees.
    height: float
        The angular height of the box in degrees.

    Returns
    -------
    box_coords: list of tuples
        The four coordinates for the box corners.

    Notes
    -----
    This calculation is only an approximation, which is not 
    self-consistent. The declination limits are calculated 
    assuming constant ra values, but the ra limits are 
    calculated at the 'lo' and 'hi' declination values. In 
    addition, this assumes the angular separation is small. 
    I found the calculation to be accurate to a few 
    percent for a box with width=height=3 degrees at a 
    declination of 80 degrees, which is fine for our purposes. 
    """
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
    box_coords: list of tuples
        The four coordinates for the box corners. This is the 
        output of the skybox function.
    butler: Butler object, optional
        If None, then a will be created in this function.
        Default is None.

    Returns
    -------
    regions: structured ndarray
        The tracts and patches within the skybox. The columns of 
        the array are 'tract' and 'patch'.
    """
    import lsst.afw.coord as afwCoord
    import lsst.afw.geom as afwGeom
    if butler is None:
        from cattools import get_butler 
        butler = get_butler()
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
