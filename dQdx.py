#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

import h5py
import os

from trackAnalyzer import *
from utils import *
from plots import *
from cuts import *

def main(infile, outdir = "None", plotfile = "None", verbose = False, plotEvents = False, plotHist = False):
    f = h5py.File(infile, 'r')
    events = f['events']
    hits = f['hits']
    tracks = f['tracks']

    tdBins = np.linspace(0, 30, 31)
    # dQdxRBins = np.linspace(0, 1.5, 31)
    dQdxBins = np.linspace(0, 50, 31)
    # dQdxRBins = np.linspace(0, 50, 51)

    tdR = np.array([])
    dQdxR = np.array([])
    td = np.array([])
    dQdx = np.array([])
    
    for track in tracks:
        if is_good_track(track):
            thisAna = trackAnalyzer(track, hits, events)
            if is_good_event(thisAna.trackEvent):

                # plot_eventDisplay(thisAna.trackEvent[0], hits, tracks)
                # plt.show()
                
                # ts, dqdx = thisAna.dQdx_ratio(np.linspace(0, 300, 16))
                # tdR = np.concatenate([tdR, ts])
                # dQdxR = np.concatenate([dQdxR, dqdx])

                ts, dqdx = thisAna.dQdx(np.linspace(0, 300, 16))
                td = np.concatenate([td, ts])
                dQdx = np.concatenate([dQdx, dqdx])

    # ratioFig = plt.figure()
    # ratioAx = ratioFig.gca()
    # plt.hist2d(tdR, dQdxR, bins = (tdBins, dQdxRBins), cmap = 'Blues')
    # plt.xlabel(r'$t - t_s$ [$\mathrm{\mu}$s]')
    # plt.ylabel(r'$dQdx$ Ratio')
    # plt.colorbar()

    if plotfile != "None":
        absFig = plt.figure()
        absAx = absFig.gca()
        plt.hist2d(td, dQdx, bins = (tdBins, dQdxBins), cmap = 'Blues')
        plt.xlabel(r'$t - t_s$ [$\mathrm{\mu}$s]')
        plt.ylabel(r'$dQdx$')
        plt.colorbar()

    import scipy.optimize as opt

    # def lt(t, a, lifetime):
    #     return a * np.exp(-t/lifetime) 
    # p0 = (30, 8)
    # def lt(t, lifetime):
    #     return np.power(2, -t/lifetime) 
    # p0 = (8)
    # def lt(t, a, b, lifetime):
    #     return a + b*np.exp(-t/lifetime)
    # p0 = (2, 1, 8)
    # def lt(t, lifetime):
    #     return np.exp(-t/lifetime)
    # p0 = (10)

    # result = opt.curve_fit(lt, tdR, dQdxR, p0 = p0)
    # print (result)
    # bf = result[0]
    # err = np.sqrt(np.diag(result[1]))

    # # plot the fitted model
    # ratioAx.plot(fineTspace, lt(fineTspace, *bf), ls = '-', color = 'red')

    # # include the fit result
    # ratioAx.text(15, 1.2,
    #              r' '.join([r'$t_\mathrm{lifetime} =',
    #                         str(round(bf[-1], 2)),
    #                         r'\pm',
    #                         str(round(err[-1], 2)),
    #                         r'$[$\mathrm{\mu}$s]']),
    #              color = 'red')

    def lt(t, a, lifetime):
        return a * np.power(2, -t/lifetime) 
    p0 = (30, 8)

    result = opt.curve_fit(lt, td, dQdx, p0 = p0)
    print ("best fit:", result)
    bf = result[0]
    err = np.sqrt(np.diag(result[1]))

    if plotfile != "None":
        # plot the fitted model
        fineTspace = np.linspace(0, 30, 1000)
        absAx.plot(fineTspace, lt(fineTspace, *bf), ls = '-', color = 'red')

        # include the fit result
        absAx.text(15, 40,
                   r' '.join([r'$t_\mathrm{lifetime} =',
                              str(round(bf[-1], 2)),
                              r'\pm',
                              str(round(err[-1], 2)),
                              r'$[$\mathrm{\mu}$s]']),
                   color = 'red')

        absFig.savefig(plotfile, dpi = 300)

    # # plot the 1-D time projection
    # fig = plt.figure()
    # plt.hist(td, bins = tdBins, histtype = 'step')
    # plt.xlabel(r'$t - t_s$ [$\mathrm{\mu}$s]')
    # plt.tight_layout()

    # # plot the 1-D dQdx ratio projection
    # fig = plt.figure()
    # plt.hist(dQdxR, bins = dQdxRBins, histtype = 'step')
    # plt.xlabel(r'$dQ/dx$ Ratio')
    # plt.tight_layout()

    # # plot the 1-D dQdx projection
    # fig = plt.figure()
    # plt.hist(dQdx, bins = dQdxBins, histtype = 'step')
    # plt.xlabel(r'$dQ/dx$ ')
    # plt.tight_layout()

    # plt.show()

    if outdir != "None":
        baseOutfileName = infile.split('/')[-1].replace('.h5', '_fitResult')
        np.save(os.path.join(outdir, baseOutfileName),
                np.array([bf[-1], err[-1]]))
        
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
                        type = str)
    # parser.add_argument('--savefit', '-s',
    #                     action = 'store_true')
    parser.add_argument('--outdir', '-o',
                        default = 'None')
    parser.add_argument('--plotfile', '-p',
                        default = 'None',
                        type = str)
    parser.add_argument('-v', '--verbose',
                        action = 'store_true')
    args = parser.parse_args()

    main(**vars(args))
