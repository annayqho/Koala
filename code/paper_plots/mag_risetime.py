""" A plot of rise time (Rest frame days)
vs observer-frame peak abs mag in g-band 
"""

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lc import get_uv_lc, get_lc
from load_lum import load_lc
from astropy.cosmology import Planck15
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

datadir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"
textsize=14


def sn2018gep(ax):
    """ Rise time and peak g-band mag for SN2018gep """
    trise = 3
    plum = -19.7

    ax.scatter(
            trise, plum, marker='o', c='k') 
    ax.text(trise/1.05, plum, "ZTF18abukavn (Ic-BL)", fontsize=textsize,
            verticalalignment='center',
            horizontalalignment='right')


def ztf(ax):
    """ Rise time and peak g-band mag for the three new ZTF sources """

    # ZTF18achdidy
    # trise = 2
    # plum = -20.2
    # ax.scatter(
    #         trise+5, plum, marker='o', c='k') 
    # ax.arrow(
    #        trise+5, plum, -5, 0, color='k', head_width=0.1, head_length=0.1)
    # ax.text(trise*1.4, plum, "ZTF18achdidy (SLSN-I)", fontsize=textsize,
    #         verticalalignment='bottom', horizontalalignment='left')

    # ZTF19aapfmki
    # trise = 4.6
    # plum = -19.5
    # ax.scatter(
    #         trise, plum, marker='o', c='k') 
    # ax.text(trise*1.08, plum*1.001, "ZTF19aapfmki (Ibn)", fontsize=textsize,
    #         verticalalignment='bottom', horizontalalignment='right')


    # ZTF18abvkwla
    trise = 1.3
    plum = -20.9
    ax.scatter(
            trise, plum, marker='o', c='k') 
    ax.text(trise*1.1, plum, "ZTF18abvkwla", fontsize=textsize,
            verticalalignment='center')


def at2018cow(ax):
    """ Rise time and peak mag """
    x = 0.4
    y = -20.6
    ax.errorbar(
            x, y, yerr=0.03, c='k', marker='o')
    ax.arrow(
           x, y, -0.1, 0, color='k', head_width=0.1, head_length=0.1)
    ax.text(x, y/1.003, "18cow", fontsize=textsize,
            horizontalalignment='center',
            verticalalignment='top')
    

def iptf16asu(ax):
    """ 
    rise time, peak Mg

    I calculated this by taking the quadratic fit to the early
    light curve in Whitesides+2017, and their estimate for the peak time.
    I calculated the t_1/2, rise using that quadratic function.
    I took their measurement of the g-band max magnitude.
    """
    # from the paper: 3.97 \pm 0.19 days in g-band
    # peak is -20.4
    trise = 1.25

    # we do resolve Mg, but Lbol is a strict lower limit
    plum = -20.4

    ax.errorbar(
            trise, plum, yerr=0.1, xerr=0.19, marker='o', c='k')
    ax.text(
            trise*1.05, plum, "iPTF16asu (Ic-BL)", fontsize=textsize,
            verticalalignment='bottom',
            horizontalalignment='left')


def rest2018(ax):
    """
    plot the data from Rest (2018)
    """
    trise = 1.4
    plum = -18.8
    ax.scatter(
            trise, plum, marker='o', c='k')
    ax.text(
            trise, plum*1.005, "15K", fontsize=textsize)


def drout(ax):
    """ using rest-frame g-band """
    x = 2.9
    y = -19.5
    ax.errorbar(
           x, y, xerr=0.1, yerr=0.08, c='k', marker='o')
    ax.text(x/1.05, y, "11qr", fontsize=textsize,
            verticalalignment='center',
            horizontalalignment='right')
    

def arcavi(ax):
    """ using rest-frame g-band """
    x = 1.12
    y = -20.7
    ax.errorbar(
           x, y, yerr=0.1, c='k', marker='o')
    ax.arrow(
           x, y, -0.1, 0, color='k', head_width=0.1, head_length=0.1)
    ax.text(x*1.05, y, "04D4ec", fontsize=textsize,
        horizontalalignment='left')

    x = 3
    y = -20.08
    ax.errorbar(
           x, y, yerr=0.1, c='k', marker='o')
    ax.text(x*1.05, y, "05D2bk", fontsize=textsize,
    verticalalignment='center', horizontalalignment='left')

    x = 4.6
    y = -20.7
    ax.errorbar(
           x, y, yerr=0.1, c='k', marker='o')
    ax.text(
            x, y*1.005, "06D1hc", 
            fontsize=textsize, horizontalalignment='right')


def ibn(ax):
    """ using rest-frame g-band """
    x = 2.3
    y = -21.26
    ax.scatter(
           x, y, c='k', marker='o')
    ax.text(x*1.05, y, "iPTF15ul (Ibn)", fontsize=textsize,
        horizontalalignment='left')

    x = 0.95
    y = -19.58
    ax.scatter(
           x, y, c='k', marker='o')
    #ax.text(x*1.05, y, "1999cq (Ibn)", fontsize=textsize,
    #    horizontalalignment='left', verticalalignment='top')
    #ax.arrow(
    #       x, y, -0.1, 0, color='k', head_width=0.1, head_length=0.1)


def dougie(ax):
    """ using rest-frame g-band, which is r-band from Vinko+2015 """
    x = 3.92
    y = -22.98
    ax.scatter(
           x, y, c='k', marker='o')
    ax.text(x/1.01, y, "Dougie", fontsize=textsize,
        horizontalalignment='right', verticalalignment='top')


def sn2011kl(ax):
    """ using rest-frame g-band, which is i-band from Greiner+2015 """
    x = 2.57
    y = -20.35
    ax.scatter(
           x, y, c='k', marker='o')
    ax.text(x*1.01, y, "SN2011kl (GRB-SN)", fontsize=textsize,
        horizontalalignment='left', verticalalignment='top')


fig,ax = plt.subplots(1,1,figsize=(6,5))
ibn(ax)
#sn2018gep(ax)
ztf(ax)
at2018cow(ax)
iptf16asu(ax)
#drout(ax)
arcavi(ax)
#dougie(ax)
sn2011kl(ax)

ax.set_ylabel("Peak Mag (Rest-frame $g$-band)", fontsize=16)
ax.set_xlim(0,5)
ax.set_ylim(-21.4, -19.9)
#ax.set_xscale('log')
ax.invert_yaxis()
ax.set_xlabel(
    r"$t_{\mathrm{1/2,rise}}$ (rest frame days)", fontsize=16)

ax.xaxis.set_tick_params(labelsize=14)
ax.yaxis.set_tick_params(labelsize=14)

fig.tight_layout()
plt.savefig("lum_rise.png", dpi=200)

#plt.show()
