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

    l_lower = 0.
    l_upper = 520. # > Detector diagonal
    l_binWidth = 10.
    l_nbin = int(np.ceil((l_upper - l_lower)/l_binWidth))
    summedCharges = np.zeros(l_nbin)
    min_axisPos = min(trackHits["axisPos"])
    for thisHit in trackHits:
        l = thisHit["axisPos"] - min_axisPos # l is in [0, length]
        l_ibin = int(np.floor(l/l_binWidth))
        summedCharges[l_ibin] += thisHit['q']
    summedCharges *= sinPolar

    if not args.use_absolute:
        norm = summedCharges[0]
        if norm == 0:
            return -1
        summedCharges /= norm

    for ibin in range(l_nbin):
        if summedCharges[ibin] > 0:
            medianL = (l_lower + l_binWidth * (ibin + 0.5))
            medianT = medianL * cosPolar / v_drift
            driftTime.append(medianT)
            charge.append(summedCharges[ibin])


    if args.use_absolute:
        bins = (np.linspace(0,200,51),
                np.logspace(1,3, 50))
        plt.hist2d(driftTime, charge,
                   #norm = LogNorm(),
                   bins = bins,
                   cmap = plt.cm.Blues)
        plt.ylabel(r'sin-corrected absolute Charge [arb.]')
    else:
        bins = (np.linspace(4,200,50), #first bin is always 1
                np.logspace(-2, 1, 50))
        plt.hist2d(driftTime, charge,
                   #norm = LogNorm(),
                   bins = bins,
                   cmap = plt.cm.Blues)
        plt.ylabel(r'Relative Charge')

    # plot a reference line
    fineTimeSpace = np.linspace(0, 200, 1000)
    if args.use_absolute:
        finePredSpace = [201*np.exp(-xi/72) for xi in fineTimeSpace]
        string = r'201exp(-t_D/72)'
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
