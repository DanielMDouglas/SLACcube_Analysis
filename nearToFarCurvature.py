#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits import mplot3d

from scipy.optimize import fmin

import h5py
import os

from wpca import WPCA

from utils import *
from plots import *

from selection import selection_dtype

def pareTrack(trackHits, trackRadialTol = 3):

    pos = np.array([trackHits['x'],
                    trackHits['y'],
                    trackHits['z']]).T
    
    CoM = np.mean(pos, axis = 0)
    
    pca = WPCA(n_components = 3)
    weight = np.stack(3*[trackHits['q']]).T
    pca.fit(pos,
            weights = weight)

    axis = pca.components_[0]
    
    axisPos = np.dot(pos-CoM,
                     axis)
    radPos = np.linalg.norm(np.cross(pos-CoM,
                                     axis), axis = -1)

    if max(radPos) > trackRadialTol:
        # print ("reducing the track by one hit!")
        # remove the most extreme point from the track
        reducedTrackHits = trackHits[~(radPos == np.max(radPos))]
        return pareTrack(reducedTrackHits)

    return pos, CoM, axis

def main(args):
    # collected all passed input files into an array
    hitData = np.empty(shape = (0), dtype = selection_dtype)
    for infile in args.infileList:
        with h5py.File(infile, 'r') as f:
            hitData = np.concatenate([hitData, f['hits']])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')

    draw_boundaries(ax)
    
    ax.set_xlabel(r'x [mm]')
    ax.set_ylabel(r'y [mm]')
    ax.set_zlabel(r'z [mm]')
    
    plt.tight_layout()
            
    trackIDs = np.unique(hitData['trackID'])
    for thisTrackID in [trackIDs[7]]:
        trackHits = hitData[hitData['trackID'] == thisTrackID]

        nearHits = trackHits[trackHits['z'] <= 50]
        farHits = trackHits[trackHits['z'] > 50]
        
        ax.scatter(farHits['x'],
                   farHits['y'],
                   farHits['z'],
                   c = "green")

        pos, CoM, axis = pareTrack(nearHits, trackRadialTol = 3)

        farPos = np.array([farHits['x'],
                           farHits['y'],
                           farHits['z'],
                           ]).T
        farRelPos = farPos - CoM
        farAxDist = np.dot(farRelPos, axis)
        farAxPos = CoM + np.outer(farAxDist, axis)
        farRadDisp = farRelPos - farAxPos

        for fAP, fP in zip(farAxPos, farPos):
            plt.plot(*np.array([fAP, fP]).T,
                     color = 'red')

        ax.scatter(*pos.T,
                   c = 'blue',
                   )

        length = 300.

        plt.plot(*np.array([CoM - axis*length/2, 
                            CoM + axis*length/2]).T,
                 color = 'black')

    plt.show()

    return

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('infileList',
                        nargs = '+',
                        type = str,
                        help = 'flattened hit data from selection.py')
    parser.add_argument('--plotfile', '-p',
                        default = "",
                        type = str,
                        help = 'optional file to which resulting lifetime measurement is saved')
    parser.add_argument('--use_absolute', '-a',
                        action = 'store_true',
                        help = 'flag to force the use of absolute charge (default: relative charge)')
    args = parser.parse_args()

    main(args)
