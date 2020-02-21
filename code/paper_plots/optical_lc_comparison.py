""" Compare the optical LC of The Koala of the light curves of
iPTF15ul, 18cow, 16X1eho, and 16asu """

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import numpy as np
from astropy.io import ascii


def iptf15ul(ax):
    dat = ascii.read("../../data/iptf15ul/iptf15ul.txt")
    mjd = dat['col2']
    filt = dat['col3']
    mag = dat['col4']
    emag = dat['col5']
    choose = np.logical_and(filt == 'g', mag != '>')


def koala(ax):
    dat = ascii.read("../../data/ZTF18abvkwla_lct.csv")
    filts = np.array([val[2:3] for val in dat['filter']])

    # Plot the g-band light curve
    choose = np.logical_and(filts=='g', dat['mag'] > 0)
    jd = dat['jdobs'][choose]
    t0 = jd[0]
    dt = (jd-t0)/(1+z)
    gmag = dat['mag'][choose] - 0.167 #extinction in SDSS g
    egmag = dat['mag_unc'][choose]
    ax.errorbar(
            dt, gmag, yerr=egmag, fmt='s', c=gcol, lw=1, ms=8)
    ax.plot(dt, gmag, c='k', lw=1)


if __name__=="__main__":
    fig,axarr = plt.subplots(2,2) 
    for ax in axarr.reshape(-1): 
        koala(ax) 
    plt.show()
