#!/usr/bin/env python 

import numpy as np
import toolbox as tb
import lsst.afw.display.ds9 as ds9
import lsst.daf.persistence
DATA_DIR = "/tigress/HSC/HSC/rerun/production-20151224/"
butler = lsst.daf.persistence.Butler(DATA_DIR)

def view(cat=None, exp=None, tract=None, patch=None, ID=None, coords=None, filter='I', frame=0, scale="zscale",\
         zoom="to fit", trans=80, draw_ells=True, maxsep=None, shape_model='shape.hsm.moments', pcolor=ds9.GREEN, ccolor=ds9.RED):

    if cat is None:
        assert (tract is not None) and (patch is not None), 'if no cat is given, must give tract and patch'
        exp = butler.get("deepCoadd_calexp", tract=tract, patch=patch, filter='HSC-'+filter, immediate=True)
        cat = butler.get("deepCoadd_meas", tract=tract, patch=patch, filter='HSC-'+filter, immediate=True)
    else:
        assert exp is not None, 'if cat is given, must also give exp'
    x0, y0 = exp.getXY0()
    settings = {'scale':scale, 'zoom': zoom, 'mask' : 'transparency %d' %(trans)}
    ds9.mtv(exp, frame=frame, settings=settings)

    #########################################################
    # maxsep is used if you only want to draw ellipses around
    # objects within a distance maxsep from a candidate
    #########################################################
    if maxsep is not None:
        if ID is not None:
            obj = cat.find(ID)
            ra0, dec0 = obj.getRa().asDegrees(), obj.getDec().asDegrees()
        else:
            ra0, dec0 = coords
        ra, dec = cat.get('coord.ra')*180./np.pi, cat.get('coord.dec')*180./np.pi
        seps = tb.angsep(ra0, dec0, ra, dec)
        cut = seps < maxsep
        cat = cat[seps < maxsep].copy(deep=True)

    if draw_ells:
        with ds9.Buffering(): 
            for i,source in enumerate(cat):
                ellcolor = pcolor if source.get('parent')==0 else ccolor
                ixx, ixy, iyy = [source.get(shape_model+_x) for _x in ['.xx', '.xy', '.yy']]
                symbol = "@:{ixx},{ixy},{iyy}".format(ixx=ixx, ixy=ixy, iyy=iyy)
                ds9.dot(symbol, source.getX()-x0, source.getY()-y0, ctype=ellcolor, frame=frame)
