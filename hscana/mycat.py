import cuts
import numpy as np

class MyCat:
    """
    This class builds a catalog of objects within a tract and patch that 
    are near a group from Yang et al. 2007. Properties such as physical 
    size and absolute magnitude are calculated assuming the redshift to
    the group. Currently, the cuts are stored as dictionaries in cuts.py.

    Initialization Parameters
    -------------------------
    group_id : int, 
        The group id from the Yang et al. 2007 catalog
    tract : int,
        HSC tract number.
    patch : string
        HSC patch. e.g., '5,7'.
    band : string, optional
        The photometric band of the observation ('G', 'R', 'I', 'Z', or 'Y').
    usewcs : bool
        If True, use the WCS to calculate the angular sizes.
    makecuts : bool
        If True, make all cuts to the catalog during the initialization

    Note: The kwargs may be used for the optional arguments to the 
          cattools.py functions.
    """
    def __init__(self, group_id, tract, patch, band='I', usewcs=False, makecuts=False, **kwargs):
        import cattools

        #######################################################
        # get group info: angular diameter distance, luminosity
        # distance, and redshift
        #######################################################

        group_info = np.genfromtxt('/home/jgreco/data/group_info.csv', delimiter=',',\
                                   dtype='i8,f8,f8,f8,f8,f8,f8,f8,i8', names=True)
        mask = group_info['group_id'] == group_id
        self.D_A = group_info['D_A'][mask][0]
        self.D_L = group_info['D_L'][mask][0]
        self.z = group_info['group_z'][mask][0]

        #######################################################
        # Get catalog and exposure for this tract, patch & band
        # if pixcuts=True, make bad pixel cuts 
        #######################################################

        butler = cattools.get_butler()
        self.exp = cattools.get_exp(tract, patch, band, butler)
        self.wcs = self.exp.getWcs() if usewcs else None
        self.cat = cattools.get_cat(tract, patch, band, butler)
        self.count_record = [] # record of number of objects
        self.count(update_record=True)

        #######################################################
        # Calculate size, mag, absolute mad, and surface 
        # brightness for all objects in the catalog.
        #######################################################

        self.angsize = cattools.get_angsize(self.cat, wcs=self.wcs, **kwargs)
        self.size = self.angsize*self.D_A*(1.0/206265.)*1.0e3 # size in kpc
        self.mag = cattools.get_mag(self.cat, self.exp.getCalib(), **kwargs)
        self.abs_mag = cattools.get_abs_mag(self.D_L, mag=self.mag)
        self.SB = cattools.get_SB(mag=self.mag, angsize=self.angsize)
        self.ra = self.cat.get('coord.ra')*180.0/np.pi
        self.dec = self.cat.get('coord.dec')*180.0/np.pi

        if makecuts:
            self.make_cuts()

    def coord(self):
        """
        Get coordinates for all objects in catalog.

        Returns
        -------
        coords : ndarry, shape = (N objects, 2)
            Coordinates (ra, dec) in degrees.
        """
        return np.dstack((self.ra, self.dec))[0]

    def count(self, update_record=False):
        """
        Count number of objects in current catalog.
        
        Parameter
        ---------
        update_record : bool, optional
            If True, update the count record, which is a list
            that saves the number of objects after each cut.

        Returns
        -------
        counts : int, only returns if update_record=False
            The number of objects in the current catalog.
        """
        num = len(self.cat)
        if update_record:
            self.count_record.append(num)
        else:
            return num

    def apply_cuts(self, cut):
        """
        Apply cuts in cut = ndarray of bools to the catalog and 
        all derived properties, and update the count record.
        """
        self.cat = self.cat[cut].copy(deep=True)
        self.size = self.angsize[cut]
        self.angsize = self.angsize[cut]
        self.mag = self.mag[cut]
        self.abs_mag = self.abs_mag[cut]
        self.SB = self.SB[cut]
        self.ra = self.ra[cut]
        self.dec = self.dec[cut]
        self.count(update_record=True)


    def make_cuts(self):
        """
        Build the cut mask and apply cuts with above method. We keep
        two records as dictionaries:

        cut_record : a record of how many objects get cut by each cut 
        nan_record : a record of how many objects get cut due to a 
            derived property (e.g., size, magnitude) begin NaN. 
        """
        from cuts import cat_cuts, phy_cuts

        # source and bad pixel cuts: the "catalog cuts"
        self.cut_record = {}
        cut  = np.ones(len(self.cat), dtype=bool)
        for col, val in cat_cuts.iteritems():
            if val is not None:
                _c = self.cat[col] == val
                self.cut_record.update({col:(~_c).sum()})
                cut &= _c
        self.apply_cuts(cut)

        # size, Mag, and SB cuts: the "physical cuts"
        self.nan_record = {}
        self.nan_record.update({'size':np.isnan(self.size).sum()})
        self.nan_record.update({'mag':np.isnan(self.mag).sum()})
        if phy_cuts['size_min'] is not None:
            self.size[np.isnan(self.size)] = -999.
            cut = self.size > phy_cuts['size_min']
            self.apply_cuts(cut)
            self.cut_record.update({'size_min':(~cut).sum()})
        if phy_cuts['SB_min'] is not None:
            self.SB[np.isnan(self.SB)] = -999.
            cut = self.SB > phy_cuts['SB_min']
            self.apply_cuts(cut)
            self.cut_record.update({'SB_min':(~cut).sum()})
        if phy_cuts['SB_max'] is not None:
            cut = self.SB < phy_cuts['SB_max']
            self.apply_cuts(cut)
            self.cut_record.update({'SB_max':(~cut).sum()})
        if phy_cuts['absmag_max'] is not None:
            cut = self.abs_mag < phy_cuts['absmag_max']
            self.apply_cuts(cut)
            self.cut_record.update({'absmag_max':(~cut).sum()})

if __name__=='__main__':
    group_id = 8480
    tract = 9349
    patch = '4,5'
    mycat = MyCat(group_id, tract, patch)
