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
        trackData = f['track'][:]

    trackHits = hits[hits['trackID'] == args.trackID]
    if( len(trackHits) == 0 ):
        return -1

    cosPolar = trackData[trackData["trackID"] == args.trackID]["cosPolar"][0]
    sinPolar = np.sqrt(1-cosPolar*cosPolar)

    driftTime = []
    charge = []

    t_lower = 0
    t_upper = 250
    t_nbin = 50
    t_binWidth = (t_upper - t_lower) / t_nbin
    # Drift time bin: [t_lower, t_upper] evenly divided into t_nbin
    # Sum up charges at same t
    # Each component is total charge at t in that bin
    summedCharges = np.zeros(t_nbin)
    for thisHit in trackHits:
        t = thisHit['z']/v_drift
        t_ibin = int(np.floor(t/t_binWidth))
        summedCharges[t_ibin] += thisHit['q']
    summedCharges *= sinPolar

    if not args.use_absolute:
        norm = summedCharges[0]
        if norm == 0:
            return -1
        summedCharges /= norm

    for ibin in range(t_nbin):
        if summedCharges[ibin] > 0:
            medianT = t_lower + t_binWidth * (ibin + 0.5)
            driftTime.append(medianT)
            charge.append(summedCharges[ibin])


    if args.use_absolute:
        bins = (np.linspace(t_lower, t_upper, t_nbin),
                np.logspace(1,3, 50))
        plt.hist2d(driftTime, charge,
                   #norm = LogNorm(),
                   bins = bins,
                   cmap = plt.cm.Blues)
        plt.ylabel(r'Absolute Charge [arb.]')
    else:
        bins = (np.linspace(t_lower, t_upper, t_nbin),
                np.logspace(-2, 1, 50))
        plt.hist2d(driftTime, charge,
                   #norm = LogNorm(),
                   bins = bins,
                   cmap = plt.cm.Blues)
        plt.ylabel(r'Relative Charge')

    # plot a reference line
    fineTimeSpace = np.linspace(t_lower, t_upper, 1000)
    if args.use_absolute:
        finePredSpace = [280*np.exp(-xi/64) for xi in fineTimeSpace]
        string = r'280exp(-t_D/64)'
    else:
        finePredSpace = [np.exp(-xi/85) for xi in fineTimeSpace]
        string = r'exp(-t_D/85)'

    plt.plot(fineTimeSpace, finePredSpace, ls = '--', color = 'red')
    
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
