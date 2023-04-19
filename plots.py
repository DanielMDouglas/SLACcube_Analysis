# some plotting utilities

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

# matplotlib.rc('font', family = 'FreeSerif', size = 16, weight = 'bold')
# matplotlib.rc('text', usetex = True)
matplotlib.rc('axes', prop_cycle = matplotlib.cycler(color = DUNEcolors))

def plot_chargeHist(q):
    chargeHistFig = plt.figure()
    chargeHistAx = chargeHistFig.gca()

    chargeHistAx.hist(q, histtype = 'step', bins = np.linspace(-50, 150, 100))
    chargeHistAx.set_xlabel(r'q')

    return chargeHistFig

def draw_hits_in_event_window_by_reference(event):
    """
    This is another way of associating hits and events
    you can use the 'ref' table to make a mask instead of t0
    """
    eventID = event['id']

    t0 = event['ts_start']

    eventMask = eventHitRefs[:,0] == eventID

    eventHits = hitData[eventMask]

    tracks = track_finder(eventHits, t0)
    
    px = eventHits['px']
    py = eventHits['py']
    ts = eventHits['ts']

    z = drift_distance(ts - t0)

    q = eventHits['q']

    X = np.array([px, py, z]).T 
    clustering = DBSCAN(eps = dbscanEps,
                        min_samples = 5).fit(X)
    # print (clustering.labels_)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')

    ax.scatter(px, py, z,
               c = clustering.labels_)
    for thisTrack in tracks:
        thisTrack.draw(ax)
        thisTrack.get_first_hit()
        print (thisTrack.is_z_fixed(10))

    return ax

def draw_boundaries(ax):
    """
    Draw the detector boundaries as a wireframe
    not needed, but pretty
    """
    boundKwargs = {'color': 'black',
                   'ls': '--'}
    
    ax.plot([detector_bounds[0][0], detector_bounds[0][1]],
            [detector_bounds[1][0], detector_bounds[1][0]],
            [detector_bounds[2][0], detector_bounds[2][0]],
            **boundKwargs)
    ax.plot([detector_bounds[0][0], detector_bounds[0][1]],
            [detector_bounds[1][1], detector_bounds[1][1]],
            [detector_bounds[2][0], detector_bounds[2][0]],
            **boundKwargs)
    ax.plot([detector_bounds[0][0], detector_bounds[0][1]],
            [detector_bounds[1][0], detector_bounds[1][0]],
            [detector_bounds[2][1], detector_bounds[2][1]],
            **boundKwargs)
    ax.plot([detector_bounds[0][0], detector_bounds[0][1]],
            [detector_bounds[1][1], detector_bounds[1][1]],
            [detector_bounds[2][1], detector_bounds[2][1]],
            **boundKwargs)

    ax.plot([detector_bounds[0][0], detector_bounds[0][0]],
            [detector_bounds[1][0], detector_bounds[1][1]],
            [detector_bounds[2][0], detector_bounds[2][0]],
            **boundKwargs)
    ax.plot([detector_bounds[0][1], detector_bounds[0][1]],
            [detector_bounds[1][0], detector_bounds[1][1]],
            [detector_bounds[2][0], detector_bounds[2][0]],
            **boundKwargs)
    ax.plot([detector_bounds[0][0], detector_bounds[0][0]],
            [detector_bounds[1][0], detector_bounds[1][1]],
            [detector_bounds[2][1], detector_bounds[2][1]],
            **boundKwargs)
    ax.plot([detector_bounds[0][1], detector_bounds[0][1]],
            [detector_bounds[1][0], detector_bounds[1][1]],
            [detector_bounds[2][1], detector_bounds[2][1]],
            **boundKwargs)

    ax.plot([detector_bounds[0][0], detector_bounds[0][0]],
            [detector_bounds[1][0], detector_bounds[1][0]],
            [detector_bounds[2][0], detector_bounds[2][1]],
            **boundKwargs)
    ax.plot([detector_bounds[0][0], detector_bounds[0][0]],
            [detector_bounds[1][1], detector_bounds[1][1]],
            [detector_bounds[2][0], detector_bounds[2][1]],
            **boundKwargs)
    ax.plot([detector_bounds[0][1], detector_bounds[0][1]],
            [detector_bounds[1][0], detector_bounds[1][0]],
            [detector_bounds[2][0], detector_bounds[2][1]],
            **boundKwargs)
    ax.plot([detector_bounds[0][1], detector_bounds[0][1]],
            [detector_bounds[1][1], detector_bounds[1][1]],
            [detector_bounds[2][0], detector_bounds[2][1]],
            **boundKwargs)

    return ax

def draw_labels(ax):
    ax.set_xlabel(r'x [mm]')
    ax.set_ylabel(r'y [mm]')
    ax.set_zlabel(r'z (drift) [mm]')

    plt.tight_layout()
