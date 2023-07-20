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

def main(args):
    # collected all passed input files into an array
    hitData = np.empty(shape = (0), dtype = selection_dtype)
    trackData = np.empty(shape = (0), dtype = track_dtype)
    for infile in args.infileList:
        with h5py.File(infile, 'r') as f:
            hitData = np.concatenate([hitData, f['hits']])
            trackData = np.concatenate([trackData, f['track']])

    driftTime = []
    absoluteCharge = []
    relativeCharge = []

    # for each track, get the initial hit and the subsequent hits
    # dQ/dz is measured from the subsequent hits using either
    # absolute or relative (to the initial hit) charge
    #trackIDs = np.unique(hitData['trackID'])
    goodTracks = trackData[(trackData["colinear"]<0.02) & (trackData["length"]>150)]
    trackIDs = goodTracks["trackID"]
    print("use", len(trackIDs), "good tracks")

    for thisTrackID in trackIDs:
        trackHits = hitData[hitData['trackID'] == thisTrackID]

        firstHit = trackHits[trackHits['z'] == 0][0]
        laterHits = trackHits[trackHits['z'] != 0]

        cosPolar = trackData[trackData["trackID"] == thisTrackID]["cosPolar"][0]
        sinPolar = np.sqrt(1-cosPolar*cosPolar)

        for thisHit in laterHits:
            
            driftTime.append(thisHit['z']/v_drift)
            relativeCharge.append(thisHit['q']/firstHit['q'])
            absoluteCharge.append(thisHit['q']*sinPolar)

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

    # fit the observation to a decay model

    if args.use_absolute:
        thisModel = expModelWithNorm((100, 100))
        thisModel.fit((driftTime, absoluteCharge))
    else:
        thisModel = expModel((100,))
        thisModel.fit((driftTime, relativeCharge))

    # plot the best fit
    fineTimeSpace = np.linspace(0, 200, 1000)
    finePredSpace = [thisModel.bf(xi) for xi in fineTimeSpace]

    plt.plot(fineTimeSpace, finePredSpace, ls = '--', color = 'red')
    plt.text(bins[0][0] + 5,
             bins[1][-1] /2.,
             thisModel.string(), color = 'red')
    plt.text(bins[0][0] + 10,
             bins[1][-1] /2.5,
             thisModel.string_errs(), color = 'red')
    
    plt.semilogy()
    plt.xlabel(r'Drift Time [$\mu$s]')
    plt.colorbar()

    plt.tight_layout()

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
    args = parser.parse_args()

    main(args)
