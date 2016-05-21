import numpy as np

cat_cuts = {'flags.pixel.bad':0,
            'flags.pixel.edge':0,
            'flags.pixel.interpolated.any':0,
            'flags.pixel.cr.any':0,
            'flags.pixel.saturated.any':0,
            'classification.extendedness':1.0,
            'parent':0,
            'flags.pixel.bright.object.any':0}

phy_cuts = {'size_min':1.0, # kpc
            'SB_min':24.0, 
            'SB_max':30.0,
            'absmag_max':-13.0}

def badpix_cut(cat):
    cut  = cat['flags.pixel.bad'] == 0
    cut &= cat['flags.pixel.edge'] == 0
    cut &= cat['flags.pixel.interpolated.any'] == 0
    cut &= cat['flags.pixel.cr.any'] == 0
    cut &= cat['flags.pixel.saturated.any'] == 0
    return cut

def source_cut(cat):
    cut  = cat['classification.extendedness'] > 0.5
    cut &= cat['parent'] == 0 
    cut &= cat['flags.pixel.bright.object.any'] == 0
    return cut

def basic_cut(cat):
    cut  = badpix_cut(cat) 
    cut &= source_cut(cat)
    return cut

def size_cut(angsize, D_A, minsize=1.0):
    """
    input parameters
    ----------------
    angsize: ndarray, angular sizes in arcsec
    D_A: float, the angular diameter distance to the group in Mpc

    optional input
    --------------
    minsize: float, minumum size for cut in kpc
    """
    size = angsize*D_A*(1.0/206265.)*1.0e3 # size in kpc
    cut = angsize > minsize
    return cut

def SB_cut(SB, SBmin=24., SBmax=30.0):
    cut = (SB>SBmin) & (SB<SBmax)
    return cut
