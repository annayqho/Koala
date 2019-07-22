""" Color evolution, luminosity evolution phase space """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from astropy.io import ascii
from astropy.cosmology import Planck15
from astropy.time import Time


def at2018gep(ax):
    distmod = Planck15.distmod(z=0.03154).value
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"
    t0 = 2458370.6473

    # For this, only use the P48 g & r forced photometry
    f = DATA_DIR + "/precursor.csv"
    prog = ascii.read(f)
    mjd = prog['mjd']
    jd = mjd + 2400000.5
    dt = jd - t0
    filt = prog['filter']
    mag = prog['mag']
    emag = prog['magerr']
    code = prog['instrument']
    det = np.logical_and.reduce((dt > 0, ~np.isnan(mag), code=='ZTF Camera'))

    # Interpolate everything onto the r-band light curve
    band = filt == 'ztfr'
    choose = np.logical_and(det, band)
    order = np.argsort(dt[choose])
    dt_r = dt[choose][order]
    mag_r = mag[choose][order]

    band = filt=='ztfg'
    choose = np.logical_and(det, band)
    order = np.argsort(dt[choose])
    dt_g = dt[choose][order]
    mag_g = mag[choose][order]

    tgrid = np.copy(dt_g)
    # keep all the points below 1 days, and then only choose 1 pt per day
    tgrid = np.hstack((tgrid[tgrid<1], np.arange(1,30,1)))
    #tgrid = np.logspace(np.log10(tgrid[0]), np.log10(tgrid[-1]), 20)

    r = np.interp(tgrid, dt_r, mag_r)
    g = np.interp(tgrid, dt_g, mag_g)

    gr = g - r
    xdata = gr
    ydata = g-distmod
    ax.plot(xdata, ydata, c='k', zorder=5, lw=2)

    markt = np.array([1/24, 1, 10, 20])
    labs = np.array(["1 hour", "1 day", "10 days", "20 days"])
    cols = np.array(['#f6d746', '#e55c30', '#84206b', '#140b34'])
    for ii,t in enumerate(markt):
        prevptx = np.interp(t-0.1, tgrid, xdata)#xdata[tgrid<t][-1]
        prevpty = np.interp(t-0.1, tgrid, ydata)#ydata[tgrid<t][-1]
        newptx = np.interp(t,tgrid,xdata)
        newpty = np.interp(t,tgrid,ydata)
        ax.annotate('', 
            xytext=(prevptx, prevpty),
            xy=(newptx, newpty),
            arrowprops=dict(color=cols[ii], width=1, headlength=10),
            label=labs[ii], zorder=10)
        ax.scatter(0,0,marker='^',c=cols[ii],s=100, label=labs[ii])
    ax.legend(loc='upper right')
        # ax.text(newptx, newpty, "$\Delta t$= %s" %labs[ii], fontsize=12,
        #         horizontalalignment='left', verticalalignment='top')
    ax.text(
            -0.3, -16, "SN2018gep", fontsize=14,
            horizontalalignment='right')


def at2018cow(ax):
    distmod = Planck15.distmod(z=0.014).value
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"
    f = DATA_DIR + "/18cow.txt"
    dat = np.loadtxt(f, dtype=str)
    jd = dat[:,0].astype(float)
    filt = dat[:,2]
    mag = dat[:,4].astype(float)
    emag = dat[:,5].astype(float)

    zp = 58285
    dt = jd-zp
    print(dt)

    tgrid = np.arange(dt[0],30,1)

    band = filt=='g'
    choose = band
    order = np.argsort(dt[choose])
    g = np.interp(tgrid, dt[choose][order], mag[choose][order])

    band = filt=='r'
    choose = band
    order = np.argsort(dt[choose])
    r = np.interp(tgrid, dt[choose][order], mag[choose][order])

    gr = g - r
    xdata = gr
    ydata = g-distmod
    ax.plot(
            xdata, ydata, c='k', ls='--', lw=2, zorder=2)
    markt = np.array([1.8, 10, 20])
    labs = np.array(["1 day", "10 days", "20 days"])
    cols = np.array(['#e55c30', '#84206b', '#140b34'])
    for ii,t in enumerate(markt):
        prevptx = np.interp(t-0.1, tgrid, xdata)#xdata[tgrid<t][-1]
        prevpty = np.interp(t-0.1, tgrid, ydata)#ydata[tgrid<t][-1]
        newptx = np.interp(t,tgrid,xdata)
        newpty = np.interp(t,tgrid,ydata)
        ax.annotate('', 
            xytext=(prevptx, prevpty),
            xy=(newptx, newpty),
            arrowprops=dict(color=cols[ii], width=1, headlength=10))
        # ax.text(newptx, newpty, "$\Delta t$= %s" %labs[ii], fontsize=12,
        #         horizontalalignment='left', verticalalignment='top')
    ax.text(-0.37, -20.4, "AT2018cow", fontsize=14)


def ksn2015k(ax):
    diff = Planck15.distmod(z=0.09).value
    gr = -0.17
    G = 20.01-diff
    ax.errorbar(
            gr, G, xerr=0.20, yerr=0.12, 
            fmt='v', c='#84206b', mfc='#84206b', ms=12, label=None)
    ax.text(gr, G, "KSN2015K", fontsize=14,
            horizontalalignment='right', verticalalignment='bottom')


if __name__=="__main__":
    fig,ax = plt.subplots(1,1,figsize=(7,5))
    cb = at2018gep(ax)
    at2018cow(ax)
    ksn2015k(ax)

    # New transient
    gr = -1*np.array([0.5, 0.5, 0.3, 0.35, 0.2, 0.15, 0])
    G = [-20.7, -20.8, -20.7, -20.6, -20.3, -20, -19.8]
    ax.errorbar(
            gr, G, 
            fmt='*', c='k', mfc='k', ms=12, label=None)

    # Formatting
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlabel("$g-r$, observer frame", fontsize=16)
    ax.set_ylabel("Absolute $g$-band mag, observer frame", fontsize=16)
    plt.xlim(-0.7, 1.0)
    plt.ylim(-15.8, -21)
    plt.legend(prop={'size':14})

    plt.tight_layout()

    plt.show()
    #plt.savefig("g_gr.eps", format="eps", dpi=1000)
