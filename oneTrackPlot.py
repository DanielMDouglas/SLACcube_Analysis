import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
matplotlib.use('TkAgg')

import h5py
from selection import selection_dtype
from consts import *


def main(args):
    with h5py.File(args.infile, 'r') as f:
        hits = f['hits'][:]

    trackHits = hits[hits['trackID'] == args.trackID]
    if( len(trackHits) == 0 ):
        return -1

    driftTime = []
    absoluteCharge = []
    relativeCharge = []

    firstHit = trackHits[trackHits['z'] == 0][0]
    laterHits = trackHits[trackHits['z'] != 0]

    for thisHit in laterHits:
        driftTime.append(thisHit['z']/v_drift)
        relativeCharge.append(thisHit['q']/firstHit['q'])
        absoluteCharge.append(thisHit['q'])

    if args.use_absolute:
        bins = (np.linspace(0, 200, 50),
                np.logspace(1,3, 50))
        plt.hist2d(driftTime, absoluteCharge,
                   #norm = LogNorm(),
                   bins = bins,
                   cmap = plt.cm.Blues)
        plt.ylabel(r'Absolute Charge [arb.]')
    else:
        bins = (np.linspace(0, 200, 50),
                np.logspace(-2, 1, 50))
        plt.hist2d(driftTime, relativeCharge,
                   #norm = LogNorm(),
                   bins = bins,
                   cmap = plt.cm.Blues)
        plt.ylabel(r'Relative Charge')

    # plot a reference line
    fineTimeSpace = np.linspace(0, 200, 1000)
    if args.use_absolute:
        finePredSpace = [90*np.exp(-xi/150) for xi in fineTimeSpace]
        string = r'90exp(-t_D/150)'
    else:
        finePredSpace = [np.exp(-xi/280) for xi in fineTimeSpace]
        string = r'exp(-t_D/280)'

    plt.plot(fineTimeSpace, finePredSpace, ls = '--', color = 'red')
    #plt.text(bins[0][0] + 5,
    #         bins[1][-1] /2.,
    #         string, color = 'red')
    
    plt.semilogy()
    plt.xlabel(r'Drift Time [$\mu$s]')
    #plt.colorbar()

    plt.tight_layout()
    plt.get_current_fig_manager().window.wm_geometry("+700+60")

    plt.show(block=False)
            
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
    parser.add_argument('--use_absolute', '-a',
                        action = 'store_true',
                        help = 'flag to force the use of absolute charge (default: relative charge)')
    args = parser.parse_args()

    main(args)
