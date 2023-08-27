#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits import mplot3d

#from scipy.optimize import fmin
from scipy.optimize import curve_fit

import h5py
import os
from wpca import WPCA

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
    def bf_tau(self):
        """
        Return tau value for best fit 
        """
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
    def bf_tau(self):
        return self.params[0]
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
    def bf_tau(self):
        return self.params[1]
    def bf_norm(self):
        return self.params[0]
    def string(self):
        return r'$' + str(round(self.params[0], 3))+ r' \exp(-t_D/' + str(round(self.params[1], 3)) + r')$' 
    def string_errs(self):
        return r'$\tau = ' + str(round(self.params[1], 3)) + r'\pm' + str(round(self.param_errs[1], 3)) +r'$'


driftTime = []
dQdxs = []
dQs = []
dxs = []
seg_colinears = []

def read_files(infileList):
    # collected all passed input files into an array
    hitData = np.empty(shape = (0), dtype = selection_dtype)
    trackData = np.empty(shape = (0), dtype = track_dtype)
    for infile in infileList:
        with h5py.File(infile, 'r') as f:
            hitData = np.concatenate([hitData, f['hits']])
            trackData = np.concatenate([trackData, f['track']])
    return hitData, trackData

def calc_dQdx(hitData, trackData):
    # for each track, get the initial hit and the subsequent hits
    # dQ/dz is measured from the subsequent hits using either
    # absolute or relative (to the initial hit) charge
    #trackIDs = np.unique(hitData['trackID'])
    goodTracks = trackData[(trackData["colinear"]<0.02) &
                           #(trackData["length"]>50) &
                           #(np.abs(trackData["cosPolar"])>0.5) &
                           (np.abs(trackData["cosPolar"]) * trackData["length"] > 250) ]

    trackIDs = goodTracks["trackID"]
    print("use", len(trackIDs), "good tracks")

    l_lower = 0.
    l_upper = 520. # > Detector diagonal
    l_binWidth = 30. #mm
    nsegment = int(np.ceil((l_upper - l_lower)/l_binWidth))
    # initialize
    driftTime.clear()
    dQdxs.clear()
    dQs.clear()
    dxs.clear()
    seg_colinears.clear()
    for thisTrackID in trackIDs:
        segment_hits = [np.empty(shape=(0), dtype = selection_dtype) for _ in range(nsegment)]
        trackHits = hitData[hitData['trackID'] == thisTrackID]

        cosPolar = trackData[trackData["trackID"] == thisTrackID]["cosPolar"][0]

        min_axisPos = min(trackHits["axisPos"])
        for thisHit in trackHits:
            l = thisHit["axisPos"] - min_axisPos # l is in [0, length]
            isegment = int(np.floor(l/l_binWidth))
            segment_hits[isegment] = np.append(segment_hits[isegment],thisHit)

        for thisSegment in segment_hits:
            if len(thisSegment) < 2:
                # Require 2 or more hits in a segment
                continue
            dQ = np.sum(thisSegment['q'])
            pos = np.array([thisSegment['x'], thisSegment['y'], thisSegment['z']]).T
            pca = WPCA()
            weight = np.stack(3*[thisSegment['q']]).T
            pca.fit(pos, weights = weight)
            seg_colinear = pca.explained_variance_[1] / pca.explained_variance_[0]
            seg_colinears.append(seg_colinear)
            if seg_colinear > 0.1:
                continue
            seg_axis = pca.components_[0]
            if seg_axis[2] < 0:
                # axis[z] should be positive
                seg_axis *= -1
            seg_axisPos = np.dot(pos, seg_axis)
            dx = np.max(seg_axisPos) - np.min(seg_axisPos)
            if np.abs(dx)<1E-9:
                continue
            dQs.append(dQ)
            dxs.append(dx)
            dQdxs.append(dQ/dx)
            zmean = np.mean(thisSegment['z'])
            tmean = zmean / v_drift
            driftTime.append(tmean)
#    return driftTime, dQdxs

def plot_lifetime():
    bins = (np.linspace(0,200,51),
            np.logspace(0,2, 51))
    plt.hist2d(driftTime, dQdxs,
               #norm = LogNorm(),
               bins = bins,
               cmap = plt.cm.Blues)
    plt.ylabel(r'dQ/dx [ADC_count/mm]')

    plt.semilogy()
    plt.xlabel(r'Drift Time [$\mu$s]')
    plt.colorbar()
    plt.tight_layout()


def detail_plots(tau):
    fig = plt.figure()

    ax_dQ = fig.add_subplot(221)
    ax_dQ.hist(dQs, bins=100)
    ax_dQ.set_xlabel("dQ [ADC_count]")
    ax_dx = fig.add_subplot(222)
    ax_dx.hist(dxs, bins=100)
    ax_dx.set_xlabel("dx [mm]")
    ax_colinear = fig.add_subplot(223)
    ax_colinear.hist(seg_colinears, bins=100)
    ax_colinear.set_xlabel("colinear parameter")

    ax_dQdx = fig.add_subplot(224)
    projected_dQdx = []
    for thisTime, thisdQdx in zip(driftTime, dQdxs):
        projected_dQdx.append(thisdQdx * np.exp(thisTime/tau))

    hist_range = (0.,100.)
    ax_dQdx.hist(dQdxs, bins=100, range=hist_range, alpha=0.5, label="Raw")
    hist, xedges, _ = ax_dQdx.hist(projected_dQdx, bins=100, range=hist_range, alpha=0.5, label="Projected")
    ax_dQdx.set_xlabel("dQ/dx [ADC_count/mm]")
    ax_dQdx.legend()
    
    plt.tight_layout()

def get_calib_const(dQdx):
    # dQdx is ADC count, pedestal subtracted, lifetime corrected, per segment length[mm]
    # Basically y-slice of fit
    dEdx = 2.2 #MeV/cm
    recomb = 0.6669
    # dQdx = dEdx * recomb * calib_const (* exp(-t/tau))
    return dQdx / recomb / dEdx *10. # ADCcount/MeV

def main(args):
    hitData, trackData = read_files(args.infileList)
    calc_dQdx(hitData, trackData)
    plot_lifetime()

    # fit the observation to a decay model
    thisModel = expModelWithNorm((100, 100))
    thisModel.fit((driftTime, dQdxs))

    # plot the best fit
    fineTimeSpace = np.linspace(0, 200, 1000)
    finePredSpace = [thisModel.bf(xi) for xi in fineTimeSpace]

    plt.plot(fineTimeSpace, finePredSpace, ls = '--', color = 'red')
    plt.text(10,20, thisModel.string(), color = 'red')
    plt.text(10,15, thisModel.string_errs(), color = 'red')

    calib_const = get_calib_const(thisModel.bf_norm())
    print("Calibration constant :",calib_const,"ADCcount/MeV")

    if args.plotfile:
        plt.savefig(args.plotfile)
    else:
        detail_plots(thisModel.bf_tau())
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
    args = parser.parse_args()

    main(args)
