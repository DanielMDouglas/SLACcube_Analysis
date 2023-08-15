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

l_binWidth = 10. #mm

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
    goodTracks = trackData[(trackData["colinear"]<0.02) &
                           #(trackData["length"]>50) &
                           #(np.abs(trackData["cosPolar"])>0.5) &
                           (np.abs(trackData["cosPolar"]) * trackData["length"] > 100) ]

    trackIDs = goodTracks["trackID"]
    print("use", len(trackIDs), "good tracks")

    l_lower = 0.
    l_upper = 520. # > Detector diagonal
    l_nbin = int(np.ceil((l_upper - l_lower)/l_binWidth))
    for thisTrackID in trackIDs:
        trackHits = hitData[hitData['trackID'] == thisTrackID]

        cosPolar = trackData[trackData["trackID"] == thisTrackID]["cosPolar"][0]
        sinPolar = np.sqrt(1-cosPolar*cosPolar)

        summedCharges = np.zeros(l_nbin)
        min_axisPos = min(trackHits["axisPos"])
        for thisHit in trackHits:
            l = thisHit["axisPos"] - min_axisPos # l is in [0, length]
            l_ibin = int(np.floor(l/l_binWidth))
            summedCharges[l_ibin] += thisHit['q']

        if use_absolute and correct_sin:
            summedCharges *= sinPolar

        if not use_absolute:
            norm = summedCharges[0]
            if norm == 0:
                print("Warning: summedCharges[0] is 0 for ID", thisTrackID)
                continue
            summedCharges /= norm

        for ibin in range(l_nbin):
            if summedCharges[ibin] > 0:
                medianL = (l_lower + l_binWidth * (ibin + 0.5))
                medianT = medianL * cosPolar / v_drift
                driftTime.append(medianT)
                charge.append(summedCharges[ibin])

    if use_absolute:
        bins = (np.linspace(0,200,51),
                np.logspace(1,3, 50))
        plt.hist2d(driftTime, charge,
                   #norm = LogNorm(),
                   bins = bins,
                   cmap = plt.cm.Blues)
        if correct_sin:
            plt.ylabel(r'sin-corrected absolute Charge [arb.]')
        else:
            plt.ylabel(r'Absolute Charge [arb.]')
    else:
        bins = (np.linspace(4,200,50), #first bin is always 1
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


def projectCharge(tau):
    fig_project = plt.figure()
    ax_project = fig_project.add_subplot(111)

    projectedCharge = []
    for thisTime, thisCharge in zip(driftTime, charge):
        if not args.use_absolute and thisCharge == 1:
            continue
        projectedCharge.append(thisCharge * np.exp(thisTime/tau))

    if args.use_absolute:
        hist_range = (0.,500.)
    else:
        hist_range = (0.,5.)
    hist, xedges, _ = ax_project.hist(projectedCharge, bins=100, range=hist_range)
    ax_project.set_xlabel("Projected Charge")

    def twoGaussian(x, *pars):
        gaus1 = pars[0] * np.exp(-1.* ((x-pars[1])/pars[2])**2)
        gaus2 = pars[3] * np.exp(-1.* ((x-pars[4])/pars[5])**2)
        return gaus1 + gaus2
    
    init_pars = [1000., hist_range[1]/3, hist_range[1]/10, 1000., hist_range[1]/3*2, hist_range[1]/10]
    hx = 0.5 * (xedges[1:]+xedges[:-1])
    result, pcov = curve_fit(twoGaussian, hx, hist, p0=init_pars)
    print("Result of two gaussian fit:")
    print("  Mean1 =",result[1]," Mean2 =",result[4])
    print("  Area1 =",result[0]*np.sqrt(np.pi)*result[2]," Area2 = ",result[3]*np.sqrt(np.pi)*result[5])
    fineX = np.linspace(*hist_range, 1000)
    fineY = [twoGaussian(xi, *result) for xi in fineX]
    ax_project.plot(fineX, fineY, ls='--', color='r')
    

def get_calib_const(ADC):
    # ADC is ADC count, pedestal subtracted, lifetime corrected, per l_binWidth[mm]
    # Basically y-slice of fit
    dEdx = 2.2 #MeV/cm
    recomb = 0.6669
    # ADC = dEdx * recomb * calib_const (* exp(-t/tau))
    return ADC / recomb / dEdx *10. # ADCcount/MeV

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

    if args.use_absolute and not args.correct_sin:
        calib_const = get_calib_const(thisModel.bf_norm())
        print("Calibration constant :",calib_const,"ADCcount/MeV")

    if args.plotfile:
        plt.savefig(args.plotfile)
    else:
        projectCharge(thisModel.bf_tau())
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
    parser.add_argument('--correct_sin', '-c',
                        action = 'store_true',
                        help = 'flag to force charge correction by sin(PolarAngle)\
                                for each track (default: no correction)')
    args = parser.parse_args()

    main(args)
