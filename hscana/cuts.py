
"""
Cuts for UDG search. The cuts are made within the MyCat class. 

cat_cuts --> catalog cuts 
phy_cuts --> physical cuts (based on parameter measurements)
"""

cat_cuts = {'flags.pixel.bad':0,
            'flags.pixel.edge':0,
            'flags.pixel.interpolated.any':0,
            'flags.pixel.cr.any':0,
            'flags.pixel.saturated.any':0,
            'classification.extendedness':1.0,
            'parent':0,
            'flags.pixel.bright.object.any':0}

phy_cuts = {'size_min':1.5, # kpc
            'SB_min':24.0, 
            'SB_max':30.0,
            'absmag_max':-13.0}
