import numpy as np

import matplotlib.pyplot as plt

import h5py

from sklearn.cluster import DBSCAN
from wpca import WPCA

from SLACplots.colors import *

import tqdm

from consts import *

selection_dtype = np.dtype([("trackID", "u4"),
                            ('q', "f4"),
                            ('x', "f4"),
                            ('y', "f4"),
                            ('z', "f4"),
                            ])

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

        self.output_arr = np.empty(len(px), dtype = selection_dtype)
        self.output_arr['q'] = q
        self.output_arr['x'] = px
        self.output_arr['y'] = py
        self.output_arr['z'] = z

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

    hitArray = np.empty(0, dtype = selection_dtype)
    for i, thisTrack in enumerate(goodTracks):
        thisTrack.output_arr['trackID'][:] = i
        hitArray = np.concatenate((hitArray, thisTrack.output_arr))

    with h5py.File(args.outfile, 'w') as of:
        of['hits'] = hitArray

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('infileList',
                        nargs = '+',
                        type = str,
                        help = 'h5flow files to parse for margin cut passing tracks')
    parser.add_argument('-o', '--outfile',
                        type = str,
                        help = 'HDF5 file to save flattened hit data from good tracks')
    args = parser.parse_args()

    main(args)
