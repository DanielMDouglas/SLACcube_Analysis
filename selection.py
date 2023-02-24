import numpy as np

import matplotlib.pyplot as plt

import h5py

from sklearn.cluster import DBSCAN
from wpca import WPCA

from SLACplots.colors import *

import tqdm

detector_bounds = [[-150, 150], [-150, 150], [0, 300]] # mm (x, y, z)

v_drift = 1.6 # mm/us (very rough estimate)
clock_interval = 0.1 # us/tick -- 10 MHz clock rate
drift_distance = detector_bounds[2][1] - detector_bounds[2][0] 

drift_window = drift_distance/(v_drift*clock_interval) # maximum drift time

drift_direction = 1 # +/- 1 depending on the direction of the drift in z

dbscanEps = 30

class track:
    def __init__(self, hits, t0):
        self.hits = hits

        px = hits['px']
        py = hits['py']
        ts = hits['ts']

        z = drift_distance(ts - t0)

        q = hits['q']

        self.pos = np.array([px, py, z]).T

        self.CoM = np.mean(self.pos, axis = 0)
        # self.axis = axis # do PCA
        pca = WPCA(n_components = 3)
        weight = np.stack(3*[q]).T
        pca.fit(self.pos, weights = weight)

        self.axis = pca.components_[0]

        self.axisPos = np.dot(self.pos-self.CoM,
                              self.axis)
        self.length = np.max(self.axisPos) - np.min(self.axisPos)
        
    def draw(self, axes):
        # print (self.CoM, self.axis)
        axes.scatter(*self.CoM, color = 'r', s = 10)
        axes.plot(*np.array([self.CoM - 0.5*self.length*self.axis,
                             self.CoM + 0.5*self.length*self.axis]).T,
                  color = 'r',
                  ls = '--')

    def get_first_hit(self):
        extremalHits = [self.pos[self.axisPos == np.min(self.axisPos)][0],
                        self.pos[self.axisPos == np.max(self.axisPos)][0]]
        # print (extremalHits)
        firstPos = min(extremalHits, key = lambda x: x[2])
        # firstZ = np.min(self.pos[:,2])
        # # print ("firstZ", firstZ)
        # firstPos = self.pos[self.pos[:,2] == firstZ][0]
        # print ("firstPos", firstPos)
        return firstPos

    def is_z_fixed(self, marginWidth):
        firstPos = self.get_first_hit()
        conditions = [any(hitPos[2] == 0 for hitPos in self.pos),
                      firstPos[0] > detector_bounds[0][0] + marginWidth,
                      firstPos[0] < detector_bounds[0][1] - marginWidth,
                      firstPos[1] > detector_bounds[1][0] + marginWidth,
                      firstPos[1] < detector_bounds[1][1] - marginWidth,
                      ]
        # print (conditions)
        # print (firstPos)
        return all(conditions)
                
        

def track_finder(hits, t0):
    px = hits['px']
    py = hits['py']
    ts = hits['ts']

    z = drift_distance(ts - t0)

    q = hits['q']

    X = np.array([px, py, z]).T 
    clustering = DBSCAN(eps = dbscanEps,
                        min_samples = 5).fit(X)
    # print (clustering.labels_)

    foundTracks = []
    for thisLabel in np.unique(clustering.labels_):
        if not thisLabel == -1:
            theseHits = hits[clustering.labels_ == thisLabel]
            foundTracks.append(track(theseHits, t0))

    return foundTracks

    # make tracks

    # for all of the found tracks, see if they are colinear
    # if they are, join them
    
def drift_distance(dt):
    """
    Estimate the z-position of a drifting electron
    """
    return detector_bounds[2][0] + drift_direction*dt*clock_interval*v_drift

def trackCut_stats_from_event(event, margin):
    eventID = event['id']

    t0 = event['ts_start']

    eventMask = eventHitRefs[:,0] == eventID

    eventHits = hitData[eventMask]

    tracks = track_finder(eventHits, t0)

    goodTracks = 0
    totalTracks = 0
    for thisTrack in tracks:
        totalTracks += 1
        if thisTrack.is_z_fixed(margin):
            goodTracks += 1

    return goodTracks, totalTracks

def good_tracks_from_event(event, hitData, eventHitRefs, margin):
    eventID = event['id']

    t0 = event['ts_start']

    eventMask = eventHitRefs[:,0] == eventID

    eventHits = hitData[eventMask]

    tracks = track_finder(eventHits, t0)

    goodTracks = []
    for thisTrack in tracks:
        if thisTrack.is_z_fixed(margin):
            goodTracks.append(thisTrack)

    return goodTracks
    
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

def main(args):

    goodTracks = []

    pbar = tqdm.tqdm(args.infileList)
    for infileName in pbar:
        pbar.set_description(infileName.split('/')[-1])
        f = h5py.File(infileName)

        eventData = f['charge']['events']['data']

        hitData = f['charge']['hits']['data']
        
        eventHitRefs = f['charge']['events']['ref']['charge']['hits']['ref']  

        for i, event in enumerate(eventData):
            goodTracks += good_tracks_from_event(event, hitData, eventHitRefs, 15)

    print ("found", len(goodTracks), "good tracks!")
    print ("saving track objects to", args.outfile)

    posData = [thisTrack.pos for thisTrack in goodTracks]
    np.savez_compressed(args.outfile, *posData)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('infileList',
                        nargs = '+',
                        type = str)
    parser.add_argument('-o', '--outfile',
                        type = str)
    args = parser.parse_args()

    main(args)

    
# marginSpace = np.linspace(0, 100, 101)
# passEff = []
# for margin in marginSpace:
#     goodTracks = 0
#     totalTracks = 0
    
#     for i, event in enumerate(eventData):
#         evGoodTracks, evTotalTracks = trackCut_stats_from_event(event, margin)

#         goodTracks += evGoodTracks
#         totalTracks += evTotalTracks

#     print (goodTracks, totalTracks, float(goodTracks/totalTracks))

#     passEff.append(float(goodTracks/totalTracks))

# plt.plot(marginSpace, passEff)
# plt.xlabel(r'Anode Veto Margin [mm]')
# plt.ylabel(r'Cut Efficiency')
# # plt.semilogy()
# plt.tight_layout()
# plt.show()

# for i, event in enumerate(eventData):
#     ax = draw_hits_in_event_window_by_reference(event)
#     draw_boundaries(ax)
#     draw_labels(ax)

#     plt.show()
