import numpy as np

def is_good_event(event):
    if event['ntracks'] > 4:
        return False
    if event['nhit'] > 1000:
        return False
    return True

def is_good_track(track):
    if track['theta'] > np.pi/2.:
        return False
    if track['length'] < 100.:
        return False
    # print (track['end'][2] - track['start'][2])
    if track['end'][2] - track['start'][2] < 40:
        return False
    return True 
