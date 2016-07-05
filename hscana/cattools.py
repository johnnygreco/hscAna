import numpy as np

def get_butler(DATA_DIR="/tigress/HSC/HSC/rerun/production-20160523/"):
    import lsst.daf.persistence
    butler = lsst.daf.persistence.Butler(DATA_DIR)
    return butler

def get_cat(tract, patch, band='I', butler=None):
    if butler is None:
        butler = get_butler()
    cat = butler.get('deepCoadd_meas', tract=tract, patch=patch, filter='HSC-'+band, immediate=True)
    return cat

def get_exp(tract, patch, band='I', butler=None):
    if butler is None:
        butler = get_butler()
    exp = butler.get('deepCoadd_calexp', tract=tract, patch=patch, filter='HSC-'+band, immediate=True)
    return exp

def get_mag(cat, calib, flux_model='cmodel.flux', **kwargs):
    calib.setThrowOnNegativeFlux(False)  # don't raise an exception when we encounter a negative or NaN flux
    mag = calib.getMagnitude(cat.get(flux_model))
    return mag

def get_absmag(D_L, cat=None, mag=None, **kwargs):
    if mag is None:
        mag = get_mag(cat, **kwargs)
    abs_mag = mag - 5.0*np.log10(D_L*1e6) + 5.0
    return abs_mag

def get_angsize(cat, wcs=None, shape_model='shape.hsm.moments', **kwargs):
    if wcs is None:
        angsize = np.power(cat.get(shape_model+'.xx')*cat.get(shape_model+'.yy')-cat.get(shape_model+'.xy')**2, 0.25)
        return angsize*0.168
    else:
        import lsst.afw.geom
        import lsst.afw.table.tableLib
        if type(cat)==lsst.afw.table.tableLib.SourceRecord:
            record = cat
            affine = wcs.linearizePixelToSky(record.get('coord'), lsst.afw.geom.arcseconds)
            shape = record.get(shape_model)
            moments = shape.transform(affine.getLinear())
            separable = lsst.afw.geom.ellipses.SeparableDistortionDeterminantRadius(moments)
            angsize = separable.getDeterminantRadius()
            return angsize
        else:
            shapekey = cat.schema.find(shape_model).key
            coordkey = cat.schema.find('coord').key
            angsize = []
            for record in cat:
                affine = wcs.linearizePixelToSky(record.get(coordkey), lsst.afw.geom.arcseconds)
                shape = record.get(shapekey)
                moments = shape.transform(affine.getLinear())
                separable = lsst.afw.geom.ellipses.SeparableDistortionDeterminantRadius(moments)
                angsize.append(separable.getDeterminantRadius())
            angsize = np.array(angsize)
            return angsize

def get_SB(cat=None, mag=None, angsize=None, return_all=False, **kwargs):
    if angsize is None:
        angsize = get_angsize(cat, **kwargs)
    if mag is None:
        mag = get_mag(cat, **kwargs)
    SB = mag + 2.5*np.log10(np.pi*angsize**2)
    return (angsize, mag, SB) if return_all else SB
