#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits import mplot3d

#from scipy.optimize import fmin
from scipy.optimize import curve_fit

import h5py
import os

from utils import *
from plots import *

from selection import selection_dtype
from selection import track_dtype

class model:
    """
    Generic model with fitting methods
    """
    def __init__(self, init_guess):
        """
        Initialize with some parameter guess, save initial guess and current bf.
        """
        self.init_params = init_guess
        self.params = init_guess
        self.param_errs = np.zeros_like(init_guess)
    def fit(self, obs):
        """
        Do a simple RMS residual minimization given an observation tuple (dep., indep.).
        """
    #    def RMS(params, x, y):
    #        residual = [yi - self(xi, params)
    #                    for xi, yi in zip(x, y)]
    #        sqResid = [pow(r, 2) for r in residual]
    #        RMSResid = pow(sum(sqResid), 0.5)
    #        return RMSResid
        
    #    result = fmin(RMS,
    #                  self.init_params,
    #                  args = (obs[0], obs[1]))
        result, pcov = curve_fit(self, obs[0], obs[1], p0=self.init_params)
        self.params = result
        self.param_errs = np.sqrt(np.diag(pcov))
    def __call__(self, x, params):
        """
        Model function.  Must implement for each model!
        """
        return x
    def bf(self, x):
        """
        Call the function with the best-fit parameters
        """
      #  return self(x, self.params)
        return 0
    def string(self):
        """
        Just a nice string with the functional form and current bf parameters
        """
        return r'$x$'
    def string_errs(self):
        """
        Just a nice string with the parameter errors
        """
        return r'$x$'

class expModel (model):
    """
    Exponential model assuming f(0) = 1
    """
    def __call__(self, x, params):
        #T, = params
        T = params
        return np.exp(-x/T)
    def bf(self, x):
        return self(x, self.params)
    def string(self):
        return r'$\exp(-t_D/' + str(round(self.params[0], 3)) + r')$' 
    def string_errs(self):
        return r'$\tau = ' + str(round(self.params[0], 3)) + r'\pm' + str(round(self.param_errs[0], 3)) +r'$'

class expModelWithNorm (model):
    """
    Exponential model including a normalization term
    """
   # def __call__(self, x, params):
   #     A, T = params
    def __call__(self, x, A, T):
        return A*np.exp(-x/T)
    def bf(self, x):
        return self(x, self.params[0], self.params[1])
    def string(self):
        return r'$' + str(round(self.params[0], 3))+ r' \exp(-t_D/' + str(round(self.params[1], 3)) + r')$' 
    def string_errs(self):
        return r'$\tau = ' + str(round(self.params[1], 3)) + r'\pm' + str(round(self.param_errs[1], 3)) +r'$'


driftTime = []
charge = []

def read_and_plot(infileList, use_absolute, correct_sin):
    # collected all passed input files into an array
    hitData = np.empty(shape = (0), dtype = selection_dtype)
    trackData = np.empty(shape = (0), dtype = track_dtype)
    for infile in infileList:
        with h5py.File(infile, 'r') as f:
            hitData = np.concatenate([hitData, f['hits']])
            trackData = np.concatenate([trackData, f['track']])

    # for each track, get the initial hit and the subsequent hits
    # dQ/dz is measured from the subsequent hits using either
    # absolute or relative (to the initial hit) charge
    #trackIDs = np.unique(hitData['trackID'])
    goodTracks = trackData[#(trackData["colinear"]<0.02) &
                           #(trackData["length"]>50) &
                           #(np.abs(trackData["cosPolar"])>0.5) &
                           (np.abs(trackData["cosPolar"]) * trackData["length"] > 100) ]

    trackIDs = goodTracks["trackID"]
    print("use", len(trackIDs), "good tracks")

    t_lower = 0
    t_upper = 250
    t_nbin = 50
    t_binWidth = (t_upper - t_lower) / t_nbin
    # Drift time bin: [t_lower, t_upper] evenly divided into t_nbin
    # Sum up charges at same t
    # Each component is total charge at t in that bin
    for thisTrackID in trackIDs:
        trackHits = hitData[hitData['trackID'] == thisTrackID]

        cosPolar = trackData[trackData["trackID"] == thisTrackID]["cosPolar"][0]
        sinPolar = np.sqrt(1-cosPolar*cosPolar)

        summedCharges = np.zeros(t_nbin)
        for thisHit in trackHits:
            t = thisHit['z']/v_drift
            t_ibin = int(np.floor(t/t_binWidth))
            summedCharges[t_ibin] += thisHit['q']

        if use_absolute and correct_sin:
            summedCharges *= sinPolar

        if not use_absolute:
            norm = summedCharges[0]
            if norm == 0:
                continue
            summedCharges /= norm

        for ibin in range(t_nbin):
            if summedCharges[ibin] > 0:
                medianT = t_lower + t_binWidth * (ibin + 0.5)
                driftTime.append(medianT)
                charge.append(summedCharges[ibin])



    if use_absolute:
        bins = (np.linspace(t_lower, t_upper, t_nbin),
                np.logspace(1,3, 50))
        plt.hist2d(driftTime, charge,
                   #norm = LogNorm(),
                   bins = bins,
                   cmap = plt.cm.Blues)
        plt.ylabel(r'Absolute Charge [arb.]')
    else:
        bins = (np.linspace(t_lower + t_binWidth, t_upper, t_nbin - 1),
                np.logspace(-2, 1, 50))
        plt.hist2d(driftTime, charge,
                   #norm = LogNorm(),
                   bins = bins,
                   cmap = plt.cm.Blues)
        plt.ylabel(r'Relative Charge')

    plt.semilogy()
    plt.xlabel(r'Drift Time [$\mu$s]')
    plt.colorbar()
    plt.tight_layout()


def main(args):
    read_and_plot(args.infileList, args.use_absolute, args.correct_sin)

    # fit the observation to a decay model
    if args.use_absolute:
        thisModel = expModelWithNorm((100, 100))
        thisModel.fit((driftTime, charge))
    else:
        thisModel = expModel((100,))
        thisModel.fit((driftTime, charge))

    # plot the best fit
    fineTimeSpace = np.linspace(0, 200, 1000)
    finePredSpace = [thisModel.bf(xi) for xi in fineTimeSpace]

    plt.plot(fineTimeSpace, finePredSpace, ls = '--', color = 'red')
    if args.use_absolute:
        text_pos = (10,20)
    else:
        text_pos = (10,5)
    plt.text(text_pos[0], text_pos[1],
             thisModel.string(), color = 'red')
    plt.text(text_pos[0], text_pos[1]/1.3,
             thisModel.string_errs(), color = 'red')

    if args.plotfile:
        plt.savefig(args.plotfile)
    else:
        plt.show()
            
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
    parser.add_argument('--correct_sin', '-s',
                        action = 'store_true',
                        help = 'flag to force charge correction by sin(PolarAngle)\
                                for each track (default: no correction)')
    args = parser.parse_args()

    main(args)
