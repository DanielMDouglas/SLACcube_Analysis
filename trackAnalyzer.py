import numpy as np

from utils import *
from consts import *

class trackAnalyzer:
    def __init__(self, track, hits, events):
        self.track = track
        self.trackHits = hits[track['hit_ref']]
        self.trackEvent = events[track['event_ref']]

        self.t0 = self.trackEvent['ts_start']
        
        
    def dQdx(self, xBins):
        startPoint = self.track['start'][:3]
        endPoint = self.track['end'][:3]

        trackDir = endPoint - startPoint
        trackDir /= mag(trackDir)

        hitPositions = np.array([self.trackHits['px'], 
                                 self.trackHits['py'],
                                 vd*clock_interval*(self.trackHits['ts'] - self.t0)]).T

        hitTs = self.trackHits['ts'] 
        hitTs -= np.min(hitTs)
        hitTs = hitTs*clock_interval
        
        trackDist = dot(hitPositions - startPoint, trackDir)
        trackDist -= trackDist[0]

        hitCharge = self.trackHits['q']
        
        sumCharge, bins = np.histogram(trackDist,
                                       weights = hitCharge,
                                       bins = xBins)
        sumTs, bins = np.histogram(trackDist,
                                   weights = hitTs,
                                   bins = xBins)
        nHits, bins = np.histogram(trackDist,
                                   bins = xBins)

        dx = 0.5*np.diff(bins)
        xMean = 0.5*(bins[1:] + bins[:-1])

        dqdx = sumCharge/dx

        meanTs = sumTs/nHits
        
        # dQdxRatio = dqdx[1:]/dqdx[0]

        # print (dqdx)
        # print (dQdxRatio.shape)

        # segTsMid = vd*clock_interval*
        lBins = startPoint + trackDir*xBins.reshape(-1, 1)
        tsBins = lBins[:,2]/vd
        tsBinCenters = 0.5*(tsBins[1:] + tsBins[:-1])
        # tsBinCenters -= tsBinCenters[0]

        # dQdxRatio = dqdx[1:]/dqdx[0]
        # tsBinCenters = tsBinCenters[1:] - tsBinCenters[0]
        dQdxRatio = dqdx
        # dQdxRatio = dqdx/dqdx[0]
        tsBinCenters = tsBinCenters - tsBinCenters[0]

        # tsBinCenters = tsBinCenters[dQdxRatio != 0]
        # dQdxRatio = dQdxRatio[dQdxRatio != 0]
        dQdxRatio = dQdxRatio[~np.isnan(meanTs)]
        meanTs = meanTs[~np.isnan(meanTs)]

        # segmentTs = tsBinCenters
        segmentTs = meanTs

        # print (dQdxRatio[0])
        return segmentTs, dQdxRatio

        # return tsBinCenters, dqdx

        # return tsBinCenters, xMean

    def dQdx_ratio(self, xBins):
        startPoint = self.track['start'][:3]
        endPoint = self.track['end'][:3]

        trackDir = endPoint - startPoint
        trackDir /= mag(trackDir)

        hitPositions = np.array([self.trackHits['px'], 
                                 self.trackHits['py'],
                                 vd*clock_interval*(self.trackHits['ts'] - self.t0)]).T

        hitTs = self.trackHits['ts']
        
        trackDist = dot(hitPositions - startPoint, trackDir)
        trackDist -= trackDist[0]

        hitCharge = self.trackHits['q']
        
        counts, bins = np.histogram(trackDist, weights = hitCharge, bins = xBins)

        dx = 0.5*np.diff(bins)
        xMean = 0.5*(bins[1:] + bins[:-1])

        dqdx = counts/dx
        
        # dQdxRatio = dqdx[1:]/dqdx[0]

        # print (dqdx)
        # print (dQdxRatio.shape)

        # segTsMid = vd*clock_interval*
        lBins = startPoint + trackDir*xBins.reshape(-1, 1)
        tsBins = lBins[:,2]/vd
        tsBinCenters = 0.5*(tsBins[1:] + tsBins[:-1])
        # tsBinCenters -= tsBinCenters[0]

        dQdxRatio = dqdx[1:]/dqdx[0]
        tsBinCenters = tsBinCenters[1:] - tsBinCenters[0]
        # dQdxRatio = dqdx
        # dQdxRatio = dqdx/dqdx[0]
        # tsBinCenters = tsBinCenters - tsBinCenters[0]

        tsBinCenters = tsBinCenters[dQdxRatio != 0]
        dQdxRatio = dQdxRatio[dQdxRatio != 0]

        # print (dQdxRatio[0])
        return tsBinCenters, dQdxRatio

        # return tsBinCenters, dqdx

        # return tsBinCenters, xMean
