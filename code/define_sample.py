""" Define a sample of transients from other surveys that satisfy
our criteria """

import numpy as np

ddir = "/Users/annaho/Dropbox/Projects/Research/Koala/data"
inputf = ddir + "/drout14.txt"
dat = np.loadtxt(inputf,delimiter=';',dtype=str)
names = np.array([val.strip() for val in dat[:,0]])
filters = np.array([val.strip() for val in dat[:,1]])
phase = np.array([val.strip() for val in dat[:,2]]).astype(float)
islim = np.array([val.strip()=='<' for val in dat[:,3]])
mag = np.array([val.strip() for val in dat[:,4]]).astype(float)

unames = np.unique(names)
ufilts = np.unique(filters)

for name in unames:
    for filt in ufilts:
        choose = np.logical_and.reduce(
                (names==name, filters==filt, islim==False))
        plt.scatter(phase[choose], mag[choose])
        plt.plot(phase[choose], mag[choose])
        indpeak = np.argmin(mag[choose])
        peak = mag[choose][indpeak]
        tpeak = phase[choose][indpeak]
        plt.axhline(y=peak+1)
        plt.gca().invert_yaxis()
        choose = np.logical_and.reduce(
                (names==name, filters==filt, islim==True))
        plt.scatter(phase[choose], mag[choose], marker='v')
