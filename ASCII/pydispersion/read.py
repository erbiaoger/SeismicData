#!/usr/bin/env python

import numpy as np
import scipy.fftpack
import matplotlib.pyplot as plt
from src import readers, writers
from src.dispersion import get_dispersion, get_fk, get_ifk


reader = getattr(readers, 'su')
writer = getattr(writers, 'su') 
def read(filename):
    traces = reader(filename)
    return traces

def read_ascii(path, NR, nt):
    from numpy import loadtxt
    from obspy.core import Stream, Stats, Trace
    dat_type = 'semd'
    comp1 = 'FXZ'
    comp2 = 'FXZ'
    stream = Stream()
    for rec_x in range(0,NR):
        file_name_in1 = path + 'P.R' + str(int(rec_x+1)) + '.' + comp1 + '.' + dat_type
        file_name_in2 = path + 'P.R' + str(int(rec_x+1)) + '.' + comp2 + '.' + dat_type
        xz1 = np.genfromtxt(file_name_in1)
        xz2 = np.genfromtxt(file_name_in2)
        deg = 0.0
        alpha = np.arctan(xz2[:nt,1]/(1.0e-40 + xz1[:nt,1])) # angle of projection
        direction = np.sign(np.cos(deg*np.pi/180.0)*xz1[:nt,1]*np.cos(alpha) + np.sin(deg*np.pi/180.0)*xz2[:nt,1]*np.cos(alpha))    
        data = direction*np.sqrt(xz1[:nt,1]**2 + xz2[:nt,1]**2)*np.cos(alpha) # scalar radial component

        stats = Stats()
        stats.filename = path + 'P.R' + str(int(rec_x+1))
        stats.starttime = xz1[0,0]
        stats.delta = xz1[1,0] - xz1[0,0]
        stats.npts = len(xz1[:nt,0])

        try:
            parts = filename.split('.')
            stats.network = parts[0]
            stats.station = parts[1]
            stats.channel = temp[2]
        except:
            pass

        stream.append(Trace(data=data[:], header=stats))

    return stream


if __name__ == '__main__':

    path = './data/guided_waves_pipe/data_L02/'
    NR = 24
    nt = 2000 # 6000
    u = read_ascii(path, NR, nt)
    dt = 3.5e-5
    dx = 10.0
    cmin = 50.
    cmax = 8000.0
    dc = 10.
    fmax = 600.0
    

    U, k, f = get_fk(u, dt, nt, dx, NR)
    print(U.shape, k.shape, f.shape)
    UU = U[:, :100]
    
    im, ax = plt.subplots(figsize=(7.0,5.0))
    plt.imshow(np.abs(UU), aspect='auto', extent=[f[0], f[100], k[0], k[-1]], cmap='jet')
    
    
    im, ax = plt.subplots(figsize=(7.0,5.0))
    traces = get_ifk(U)
    plt.imshow(np.abs(traces), aspect='auto', cmap='jet')
    
    im, ax = plt.subplots(figsize=(7.0,5.0))
    plt.imshow(u, aspect='auto', cmap='jet')
