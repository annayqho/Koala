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
textsize=12


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

    # ZTF18abvkwla
    trise = 1.83
    plum = -20.59
    ax.errorbar(
            trise, plum, xerr=0.05, yerr=0.07, marker='*', c='k', ms=20) 
    ax.text(trise, plum*1.003, "ZTF18abvkwla", fontsize=textsize+2,
            verticalalignment='bottom', weight='bold', horizontalalignment='center')


def at2018cow(ax):
    """ Rise time and peak mag """
    x = 1.43
    y = -20.87
    ax.errorbar(
            x, y, xerr=0.08, yerr=0.05, c='k', marker='o')
    ax.text(x/1.03, y, "18cow", fontsize=textsize,
            horizontalalignment='right',
            verticalalignment='top')
    

def iptf16asu(ax):
    """ 
    rise time, peak Mg

    I calculated this by taking the quadratic fit to the early
    light curve in Whitesides+2017, and their estimate for the peak time.
    I calculated the t_1/2, rise using that quadratic function.
    I took their measurement of the g-band max magnitude.
    """
    trise = 2.73
    plum = -20.2

    ax.errorbar(
            trise, plum, yerr=0.1, xerr=0.19, marker='o', c='k')
    ax.text(
            trise/1.05, plum, "iPTF16asu (Ic-BL)", fontsize=textsize,
            verticalalignment='bottom',
            horizontalalignment='right')


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
    # 04D4ec
    x = 3.81
    y = -20.26
    ax.errorbar(
           x, y, yerr=0.03, c='k', marker='o')
    ax.arrow(
           x, y, -0.1, 0, color='k', head_width=0.1, head_length=0.1)
    ax.text(x/1.06, y, "04D4ec", fontsize=textsize,
        horizontalalignment='right', verticalalignment='center')

    # 05D2bk
    x = 2.90
    y = -20.39
    ax.errorbar(
           x, y, xerr=0.06, yerr=0.02, c='k', marker='o')
    ax.text(x/1.01, y*1.001, "05D2bk", fontsize=textsize,
    verticalalignment='bottom', horizontalalignment='left')

    # 06D1hc
    x = 4.59
    y = -20.28
    ax.errorbar(
           x, y, xerr=0.06, yerr=0.03, c='k', marker='o')
    ax.text(
            x, y/1.001, "06D1hc", 
            fontsize=textsize, horizontalalignment='center',
            verticalalignment='top')


def ibn(ax):
    """ using rest-frame g-band """
    x = 1.5
    y = -21.2
    ax.errorbar(
           x, y, xerr=0.05, yerr=0.3, c='k', marker='o')
    ax.text(x*1.05, y, "iPTF15ul (Ibn?)", fontsize=textsize,
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
    y = -23.03
    ax.errorbar(
           x, y, xerr=0.14, yerr=0.13, c='k', marker='o')
    ax.text(x/1.01, y, "Dougie", fontsize=textsize,
        horizontalalignment='right', verticalalignment='top')


def sn2011kl(ax):
    """ using rest-frame g-band, which is i-band from Greiner+2015
    and Kann+2019 """
    x = 4.97
    y = -20.31
    ax.errorbar(
           x, y, xerr=1.20, yerr=0.13, c='k', marker='o')
    ax.text(x/1.01, y*1.001, "SN2011kl (GRB-SN)", fontsize=textsize,
        horizontalalignment='right', verticalalignment='bottom')


def des16X1eho(ax):
    """ using rest-frame i-band LC from Miika """
    x = 1.9
    y = -20.39
    ax.errorbar(
           x, y, xerr=0.6, yerr=0.09, c='k', marker='o')
    ax.text(x*1.01, y*1.001, "16X1eho", fontsize=textsize,
        horizontalalignment='left', verticalalignment='bottom')


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
des16X1eho(ax)

ax.set_ylabel("Peak Mag (Rest-frame $g$-band)", fontsize=16)
#ax.set_xlim(1,5)
ax.set_ylim(-21.6, -20)
#ax.set_xscale('log')
ax.invert_yaxis()
ax.set_xlabel(
    r"$t_{\mathrm{1/2,rise}}$ (rest frame days)", fontsize=16)

ax.xaxis.set_tick_params(labelsize=14)
ax.yaxis.set_tick_params(labelsize=14)

fig.tight_layout()
plt.savefig("lum_rise.eps", dpi=300)

#plt.show()
