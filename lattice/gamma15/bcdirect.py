# -*- coding: utf-8 -*-
"""Model-free Bost-Charles evaluator: the direct double sum
   BC = (1/N^2) sum_{j,k} log|(f_j-f_k)/(z_j-z_k)|,  diagonal log|f'(z_j)|,
valid because int int_{T^2} log|z-w| dmu dmu = 0.  Slow, and needs a shrink < 1 so
that the boundary curve does not traverse the zero-width slits twice."""
import numpy as np, contour as C

def BC_direct(reg, N=2048, K=2500, s=1.0):
    zs = np.exp(2j*np.pi*np.arange(N)/N)
    f = s*np.array(C.h(np.array(reg.fwd(zs)), K))
    c = np.fft.fft(f)/N
    k = np.fft.fftfreq(N, d=1.0/N)
    fp = (np.fft.ifft(c*(1j*k))*N)/(1j*zs)
    tot = 0.0
    for j in range(N):
        with np.errstate(divide='ignore', invalid='ignore'):
            q = np.abs((f-f[j])/(zs-zs[j]))
        q[j] = abs(fp[j])
        tot += np.sum(np.log(q))
    return tot/N**2
