#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np

import h5py

from plots import *

def main(args):
    with h5py.File(args.infile, 'r') as f:
        hits = f['hits'][:]

    trackHits = hits[hits['trackID'] == args.trackID]
        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')

    draw_boundaries(ax)

    ax.scatter(trackHits['x'],
               trackHits['y'],
               trackHits['z'],
               c = trackHits['q'])
    
    ax.set_xlabel(r'x [mm]')
    ax.set_ylabel(r'y [mm]')
    ax.set_zlabel(r'z [mm]')
    
    plt.tight_layout()

    if args.plotfile:
        plt.savefig(args.plotfile)
    else:
        plt.show()
        
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
                        type = str,
                        help = 'flattened hit data from selection.py')
    parser.add_argument('-t', '--trackID',
                        type = int,
                        default = 0,
                        help = 'trackID to plot')
    parser.add_argument('--plotfile', '-p',
                        default = "",
                        type = str,
                        help = 'optional file to which resulting track display is saved')
    args = parser.parse_args()

    main(args)

