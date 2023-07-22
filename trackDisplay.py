#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('TkAgg')

import h5py

from plots import *

def main(args):
    with h5py.File(args.infile, 'r') as f:
        hits = f['hits'][:]
        if args.verbose:
            tracks = f['track'][:]

    trackHits = hits[hits['trackID'] == args.trackID]
    if( len(trackHits) == 0 ):
        return -1
        
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
    ax.view_init(elev=5)
    
    plt.tight_layout()

    if args.verbose:
        thisTrack = tracks[tracks["trackID"] == args.trackID]
        print("trackID:", thisTrack["trackID"][0],
            ", colinear", round(thisTrack["colinear"][0],4),
            ", length", round(thisTrack["length"][0],1),
            ", cosPolar", round(thisTrack["cosPolar"][0],3),
            ", passCathode", thisTrack["passCathode"][0])

    if args.plotfile:
        plt.savefig(args.plotfile)
    else:
        if __name__ == '__main__':
            plt.show()
        else:
            plt.get_current_fig_manager().window.wm_geometry("+20+60")
            plt.show(block=False)

    return 0
        
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
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help = 'print track information')
    parser.add_argument('--plotfile', '-p',
                        default = "",
                        type = str,
                        help = 'optional file to which resulting track display is saved')
    args = parser.parse_args()

    main(args)

