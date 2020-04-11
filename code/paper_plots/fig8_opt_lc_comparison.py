""" Compare the optical LC of The Koala of the light curves of
iPTF15ul, 18cow, 16X1eho, and 16asu """

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import numpy as np
import pandas as pd
from astropy.io import ascii
from astropy.cosmology import Planck15


def iptf15ul(ax,c):
    z = 0.066
    dat = ascii.read("../../data/iptf15ul/iptf15ul_corr.txt")
    mjd = dat['col2']
    filt = dat['col3']
    dm = Planck15.distmod(z=z).value
    mag = dat['col4']
    emag = dat['col5']
    choose = np.logical_and(filt == 'g', mag != '>')
    dt = (mjd[choose]-mjd[choose][0])/(1+z)
    dm = -21.1-min(mag[choose].astype(float))
    ax.plot(
        dt-1.2, mag[choose].astype(float)+dm+2.5*np.log10(1+z), c=c)


def at2018cow(ax,c):
    z = 0.0141
    dat = pd.read_fwf("../../data/at2018cow_photometry_table.dat")
    mjd = dat['MJD']
    filt = dat['Filt']
    mag = dat['ABMag']
    emag = dat['Emag']
    choose = filt == 'g'

    xall = mjd[choose][1:].values
    yall = mag[choose][1:].values.astype(float)-0.287
    yerrall = emag[choose][1:].values.astype(float)

    # Add the i-band detection prior to peak
    xall = np.insert(xall, 0, 58286.1950)
    yall = np.insert(yall, 0, 14.32-0.147-0.7)
    yerrall = np.insert(yerrall, 0, 0.03)

    # Add the o-band detection prior to peak
    xall = np.insert(xall, 0, 58285.44)
    yall = np.insert(yall, 0, 14.70-0.198-0.4)
    yerrall = np.insert(yerrall, 0, 0.10)

    dt = xall-xall[0]
    Yall = yall-Planck15.distmod(z=0.0141).value
    ax.plot(dt, Yall, c=c)


def des(ax,c):
    z = 0.76
    dat = pd.read_table("../../data/DES16X1eho_i.dat")
    xall = (dat['# t'][1:-1].values.astype(float))/(1+z)
    yall= dat['M'][1:-1].values.astype(float)
    kcorr = 2.5*np.log10(1+z)
    yerrall = dat['M_e'][1:-1].values.astype(float)
    ax.plot(xall+2, yall+kcorr, c=c)


def iptf16asu(ax,c):
    z = 0.187
    dat = pd.read_fwf("../../data/iptf16asu.txt")
    filt = dat['filt']
    # choose = filt == 'r'
    # dt = dat['dt'][choose] / (1+z)
    # mag = dat['mag'][choose]
    # Mag = mag[0:10].astype(float) - Planck15.distmod(z=z).value + 2.5*np.log10(1+z)
    # ax.plot(dt[0:10], Mag, c=c)

    choose = filt == 'g'
    dt = dat['dt'][choose] / (1+z)
    mag = dat['mag'][choose][16:30].astype(float).values
    Mag = mag - Planck15.distmod(z=z).value + 2.5*np.log10(1+z)
    ax.plot(dt[16:30]+3, Mag, c=c)


def koala(ax):
    z = 0.2714
    dat = ascii.read("../../data/ZTF18abvkwla_lct.csv")
    filts = np.array([val[2:3] for val in dat['filter']])

    # Plot the r-band light curve (rest-frame g)
    choose = np.logical_and(filts=='r', dat['mag'] > 0)
    jd = dat['jdobs'][choose]
    t0 = jd[0]
    dt = (jd-t0)/(1+z)
    gmag = dat['mag'][choose] - 0.115 #extinction in SDSS r
    egmag = dat['mag_unc'][choose]
    Gmag = gmag-Planck15.distmod(z=z).value+2.5*np.log10(1+z)
    ax.errorbar(
            dt, Gmag, yerr=egmag, fmt='s', c='k', lw=1, ms=5, zorder=5)
    ax.plot(dt, Gmag, c='k', lw=3, zorder=5)
    choose = np.logical_and(filts=='g', dat['mag'] == 0)
    jd = dat['jdobs'][choose]
    dt = (jd-t0)/(1+z)
    lims = dat['limmag'][choose] - 0.167 # extinction in SDSS g
    Lims = lims - Planck15.distmod(z=z).value + 2.5*np.log10(1+z)
    ax.scatter(
            dt, Lims, edgecolor='k', facecolor='white', marker='s', label=None,
            zorder=5)
    for ii,dt_val in enumerate(dt):
        ax.arrow(dt_val, Lims[ii], 0, +0.1, color='k',
                length_includes_head=True,
                head_width=0.2, head_length=0.05, label=None, zorder=5)


def sn2011kl(ax, c):
    z = 0.677
    dat = ascii.read("../../data/2011kl.txt")
    xall = dat['col1']
    yall = dat['col2']
    yerrall = dat['col3']
    ax.plot(xall-12, yall+2.5*np.log10(1+z), c=c)


def snls05d2bk(ax, c):
    z = 0.699
    mw_ext = 0.029

    dat = pd.read_fwf("../../data/SNLS05D2bk.txt")
    jd = dat['JD']
    filt = dat['F']
    mag = dat['mag']
    emag = dat['emag']

    choose = filt == 'i'
    xall = jd[choose][0:-2].values.astype(float)
    yall = mag[choose][0:-2].values.astype(float)
    yerrall = emag[choose][0:-2].values.astype(float)
    ax.plot(
            (xall-xall[0]-4)/(1+z), 
            yall-mw_ext-Planck15.distmod(z=z).value+2.5*np.log10(1+z), 
            c=c)
    


if __name__=="__main__":
    fig,axarr = plt.subplots(3,2, sharex=True, sharey=True, figsize=(8,8)) 
    for ax in axarr.reshape(-1): 
        koala(ax) 
        iptf15ul(ax, 'lightgrey')
        at2018cow(ax, 'lightgrey')
        des(ax, 'lightgrey')
        iptf16asu(ax, 'lightgrey')
        sn2011kl(ax, 'lightgrey')
        snls05d2bk(ax, 'lightgrey')
    ax = axarr[0,0]
    iptf15ul(ax, '#d95f0e')
    ax.text(0.95,0.95,'iPTF15ul', fontsize=12,
        horizontalalignment='right', verticalalignment='top',
        transform=ax.transAxes)
    ax = axarr[0,1]
    at2018cow(ax, '#d95f0e')
    ax.text(0.95,0.95,'AT2018cow', fontsize=12,
        horizontalalignment='right', verticalalignment='top',
        transform=ax.transAxes)
    ax = axarr[1,0]
    des(ax, '#d95f0e')
    ax.text(0.95,0.95,'DES16X1eho', fontsize=12,
        horizontalalignment='right', verticalalignment='top',
        transform=ax.transAxes)
    ax = axarr[1,1]
    iptf16asu(ax, '#d95f0e')
    ax.text(0.95,0.95,'iPTF16asu', fontsize=12,
        horizontalalignment='right', verticalalignment='top',
        transform=ax.transAxes)
    ax = axarr[2,0]
    sn2011kl(ax, "#d95f0e")
    ax.text(0.95,0.95,'SN2011kl', fontsize=12,
        horizontalalignment='right', verticalalignment='top',
        transform=ax.transAxes)
    ax =axarr[2,1]
    snls05d2bk(ax, "#d95f0e")
    ax.text(0.95,0.95,'SNLS05D2bk', fontsize=12,
        horizontalalignment='right', verticalalignment='top',
        transform=ax.transAxes)

    # Labeling
    fig.text(
            0.5, 0.04, r"$\Delta t$ (rest-frame days)", 
            ha='center', fontsize=16)
    fig.text(
            0.04, 0.5, r'Absolute Mag',
            fontsize=16, rotation='vertical', horizontalalignment='center',
            verticalalignment='center')

    # Formatting
    axarr[0,0].set_xlim(-1.3,8)
    axarr[0,0].set_ylim(-21.5,-18.5)
    axarr[0,0].invert_yaxis()
    #fig.subplots_adjust(wspace=0.05, hspace=0.2, left=0.15, bottom=0.15)
    for ax in axarr.reshape(-1):
        ax.yaxis.set_tick_params(labelsize=14)
        ax.xaxis.set_tick_params(labelsize=14)

    #plt.show()
    plt.savefig(
        "optical_lc_comparison.png", dpi=500, 
        bbox_inches='tight', pad_inches=0.1)
