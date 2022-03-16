import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

from consts import *

DUNEblue = '#7FAED5'
DUNElightOrange = '#F19E54'
DUNEdarkOrange = '#F0652B'
DUNEgreen = '#8ACA6F'
DUNEgray = '#626466'
DUNEyellow = '#FBD03F'
DUNEpurple = '#5B3069'
DUNElightPurple = '#8C6E96'
DUNEcyan = '#42C2A8'
DUNEpink = '#F7AEC2'

DUNEcolors = [DUNEblue,
              DUNElightOrange,
              DUNEgreen,
              DUNEdarkOrange,
              DUNEyellow,
              DUNEpink,
              DUNEpurple,
              DUNEcyan,
              DUNEgray,
              DUNElightPurple,
]

# LaTeXflavor = {"numu": r'$\nu_\mu$',
#                "numubar": r'$\bar{\nu}_\mu$',
#                "nue": r'$\nu_e$',
#                "nuebar": r'$\bar{\nu}_e$',
#                "nutau": r'$\nu_\tau$',
#                "nutaubar": r'$\bar{\nu}_\tau$'}

matplotlib.rc('font', family = 'FreeSerif', size = 16, weight = 'bold')
# matplotlib.rc('text', usetex = True)
matplotlib.rc('axes', prop_cycle = matplotlib.cycler(color = DUNEcolors))

def plot_chargeHist(q):
    chargeHistFig = plt.figure()
    chargeHistAx = chargeHistFig.gca()

    chargeHistAx.hist(q, histtype = 'step', bins = np.linspace(-50, 150, 100))
    chargeHistAx.set_xlabel(r'q')

    return chargeHistFig

def plot_eventDisplay(event, hits, tracks):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')

    hitRef = event['hit_ref']
    eventHits = hits[hitRef]
    
    t0 = event['ts_start']
        
    x = eventHits['px']
    y = eventHits['py']
    t = clock_interval*(eventHits['ts'] - t0)
    
    q = eventHits['q']
    
    ax.scatter(x, y, t, c = q)
    ax.set_xlabel(r'x [mm]')
    ax.set_ylabel(r'y [mm]')
    ax.set_zlabel(r't [$\mu$s]')

    if event['ntracks']:
        trackRef = event['track_ref']
        eventTracks = tracks[trackRef]
            
        for track in eventTracks:
            start4vec = track['start']
            end4vec = track['end']
            
            startPoint = [start4vec[0],
                          start4vec[1],
                          clock_interval*start4vec[3]]
            endPoint = [end4vec[0],
                        end4vec[1],
                        clock_interval*end4vec[3]]

            ax.plot(*zip(startPoint[:3], endPoint[:3]),
                    ls = '--', c = 'r')

