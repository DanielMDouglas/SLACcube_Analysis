import numpy as np

def dot(a, b):
    return np.sum(a*b, axis = -1)

def mag(a):
    return np.sqrt(np.dot(a, a))
