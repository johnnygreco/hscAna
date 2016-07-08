"""
A collection of functions for getting 
and interacting with HSC images.
"""

import numpy as np
import pipeTools as pT

def write_fits(tract, patch, band='I', butler=None, save=True, outdir='/home/jgreco/projects/hscAna/output/'):
    if butler is None:
        butler = pT.get_butler()
    print 'getting '
    calexp = pT.get_exp(tract, patch, band, butler)
    outfile = outdir+'deepCoadd_HSC-'+band+'_'+str(tract)+'_'+patch+'.fits'
    calexp.writeFits(outfile)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Get and save image given its tract and patch')
    parser.add_argument('tract', type=int, help='tract of observation')
    parser.add_argument('patch', type=str, help='patch of observation')
    parser.add_argument('-b', '--band', help='observation band', default='I')
    parser.add_argument('-o', '--outdir', help='output directory', 
                        default='/home/jgreco/projects/hscAna/output/')
    args = parser.parse_args()
    write_fits(args.tract, args.patch, band=args.band, outdir=args.outdir)
