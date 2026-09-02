# -*- coding: utf-8 -*-
import numpy as np
def winding(curve, pt):
    """Winding number of a closed polygonal curve (array, first point != last) about pt."""
    v = np.asarray(curve) - pt
    a = np.angle(v)
    d = np.diff(np.concatenate([a, a[:1]]))
    d = (d + np.pi) % (2*np.pi) - np.pi
    return d.sum()/(2*np.pi)
