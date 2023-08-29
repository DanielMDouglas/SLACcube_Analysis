import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import h5py

from sklearn.cluster import DBSCAN
from wpca import WPCA
import tqdm

from consts import *

selection_dtype = np.dtype([("trackID", "u4"),
                            ('q', "f4"),
                            ('x', "f4"),
                            ('y', "f4"),
                            ('z', "f4"),
                            ("axisPos", "f4")
                            ])
track_dtype = np.dtype([("trackID", "u4"),
                        ("totalCharge", "f4"),
                        ("colinear", "f4"),
                        ("length", "f4"),
                        ("cosPolar", "f4"),
                        ("unixtime", "u8")
                        ])

dbscanEps = 30
margin = 15

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
        if self.axis[2] < 0:
            # axis[z] should be positive
            self.axis *= -1

        self.axisPos = np.dot(self.pos-self.CoM,
                              self.axis)
        self.output_arr["axisPos"] = self.axisPos

        self.length = np.max(self.axisPos) - np.min(self.axisPos)
        self.colinear = pca.explained_variance_[1] / pca.explained_variance_[0]
        z_axis = np.array([0, 0, 1])
        self.cosPolar = np.dot(self.axis, z_axis)
        self.totalCharge = np.sum(q)
        self.unixtime = 0 #set later
        
    def draw(self, axes):
        axes.scatter(*self.CoM, color = 'r', s = 10)
        axes.plot(*np.array([self.CoM - 0.5*self.length*self.axis,
                             self.CoM + 0.5*self.length*self.axis]).T,
                  color = 'r',
                  ls = '--')

    def get_first_hit(self):
        extremalHits = [self.pos[self.axisPos == np.min(self.axisPos)][0],
                        self.pos[self.axisPos == np.max(self.axisPos)][0]]
        firstPos = min(extremalHits, key = lambda x: x[2])
        return firstPos

    def is_z_fixed(self, marginWidth):
        firstPos = self.get_first_hit()
        conditions = [any(hitPos[2] == 0 for hitPos in self.pos),
                      firstPos[0] > detector_bounds[0][0] + marginWidth,
                      firstPos[0] < detector_bounds[0][1] - marginWidth,
                      firstPos[1] > detector_bounds[1][0] + marginWidth,
                      firstPos[1] < detector_bounds[1][1] - marginWidth,
                      ]
        return all(conditions)

    def is_good_track(self):
        is_colinear = self.colinear < 0.02
        is_zlong = self.length * self.cosPolar > 100
        return is_colinear and is_zlong

def track_finder(hits, t0):
    px = hits['px']
    py = hits['py']
    ts = hits['ts']

    z = drift_distance(ts - t0)

    q = hits['q']

    X = np.array([px, py, z]).T 
    clustering = DBSCAN(eps = dbscanEps,
                        min_samples = 5).fit(X)

    foundTracks = []
    for thisLabel in np.unique(clustering.labels_):
        if not thisLabel == -1:
            theseHits = hits[clustering.labels_ == thisLabel]
            foundTracks.append(track(theseHits, t0))

    return foundTracks

    
def drift_distance(dt):
    """
    Estimate the z-position of a drifting electron
    """
    return detector_bounds[2][0] + drift_direction*dt*clock_interval*v_drift

def good_tracks_from_event(event, hitData, eventHitRefs, margin):
    eventID = event['id']

    t0 = event['ts_start']

    eventMask = eventHitRefs[:,0] == eventID

    eventHits = hitData[eventMask]

    tracks = track_finder(eventHits, t0)

    goodTracks = []
    for thisTrack in tracks:
        if thisTrack.is_z_fixed(margin) and thisTrack.is_good_track():
            thisTrack.unixtime = event['unix_ts']
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
            goodTracks += good_tracks_from_event(event, hitData, eventHitRefs, margin)

    print ("found", len(goodTracks), "good tracks!")
    print ("saving track objects to", args.outfile)

    hitArray = np.empty(0, dtype = selection_dtype)
    trackArray = np.empty(0, dtype = track_dtype)
    for i, thisTrack in enumerate(goodTracks):
        thisTrack.output_arr['trackID'][:] = i
        hitArray = np.concatenate((hitArray, thisTrack.output_arr))
        trackInfo = np.array((i, thisTrack.totalCharge,
                                 thisTrack.colinear, 
                                 thisTrack.length, 
                                 thisTrack.cosPolar,
                                 thisTrack.unixtime),
                             dtype=track_dtype)
        trackArray = np.append(trackArray, trackInfo)

    # Plot track features
    fig, axes = plt.subplots(nrows=3, ncols=3)
    axes[0,0].hist(trackArray["colinear"], bins=100)
    axes[0,0].set_xlabel("colinear")
    axes[0,0].semilogy()
    axes[1,1].hist(trackArray["length"], bins=100)
    axes[1,1].set_xlabel("length")
    axes[2,2].hist(trackArray["cosPolar"], bins=100)
    axes[2,2].set_xlabel("cosPolar")
    axes[0,1].hist2d(trackArray["length"], trackArray["colinear"], cmap=plt.cm.Blues, bins=[20,50], norm = LogNorm())
    axes[0,2].hist2d(trackArray["cosPolar"], trackArray["colinear"], cmap=plt.cm.Blues, bins=[20,50], norm = LogNorm())
    axes[1,2].hist2d(trackArray["cosPolar"], trackArray["length"], cmap=plt.cm.Blues, bins=20)
    axes[1,0].axis("off")
    axes[2,0].axis("off")
    axes[2,1].axis("off")
    if args.plotfile:
        plt.savefig(args.plotfile)
    else:
        plt.show()

    with h5py.File(args.outfile, 'w') as of:
        of['hits'] = hitArray
        of['track'] = trackArray

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('infileList',
                        nargs = '+',
                        type = str,
                        help = 'h5flow files to parse for margin cut passing tracks')
    parser.add_argument('-o', '--outfile',
                        required = True,
                        type = str,
                        help = 'HDF5 file to save flattened hit data from good tracks')
    parser.add_argument('--plotfile', '-p',
                        default = "",
                        type = str,
                        help = 'optional file to which distributions of feature variables are saved')
    args = parser.parse_args()

    main(args)
