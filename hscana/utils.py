def angsep(ra1, dec1, ra2, dec2, sepunits='arcsec'):
    """
    Return angular separation btwn points

    Inputs
    ------
    ra1:  float or ndarray, right ascension(s) of first point (in degrees)
    dec1: float or ndarray, declination(s) of first point (in degrees)
    ra2:  float or ndarray, right ascension(s) of second point (in degrees)
    dec2: float or ndarray, declination(s) of second point (in degrees)

    Optional input
    --------------
    sepunits: string, angular unit of the returned separation
              (radian, degree, or arcsec; default = arcsec)

    Returns:
    sep: same type as input angles, angular separation between points
    """
    import numpy as np
    deg2rad = np.pi/180.0

    ra1 = ra1*deg2rad
    dec1 = dec1*deg2rad
    ra2 = ra2*deg2rad
    dec2 = dec2*deg2rad

    sin_dRA = np.sin(ra2 - ra1)
    cos_dRA = np.cos(ra2 - ra1)
    sin_dec1 = np.sin(dec1)
    sin_dec2 = np.sin(dec2)
    cos_dec1 = np.cos(dec1)
    cos_dec2 = np.cos(dec2)

    num1 = cos_dec2*sin_dRA
    num2 = cos_dec1*sin_dec2 - sin_dec1*cos_dec2*cos_dRA
    denom = sin_dec1*sin_dec2 + cos_dec1*cos_dec2*cos_dRA
    sep = np.arctan2(np.sqrt(num1*num1 + num2*num2), denom)

    conversion = {'radian':1.0, 'arcsec':206264.806, 'degree':180./np.pi}
    sep *= conversion[sepunits]

    return sep

def pickle_data(data, filename):
    """
    Pickle data

    input
    -----
    data : any data construct. the data
        to be pickled.
    filename : string. name of file to 
        output with full path name.
    """
    import cPickle as pickle
    pkl_file = open(filename, 'wb')
    pickle.dump(data, pkl_file)
    pkl_file.close()

def read_pickled_data(filename):
    """
    Read the pickled data

    input
    -----
    filename : string. name of file that 
        contains pickled data

    output
    ------
    the unpickled data
    """
    import cPickle as pickle
    pkl_file = open(filename, 'rb')
    data = pickle.load(pkl_file)
    pkl_file.close()
    return data
