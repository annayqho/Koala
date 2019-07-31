""" Plot absolute magnitude vs. rise time

I'm interested in knowing what peaks between -18 and -21 mag
with rise times between 1 day and 10 days
"""


import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
from astropy.cosmology import Planck15


fig, ax = plt.subplots(1,1, figsize=(6,5))

losx = 0.1
losy = -0.06


def plot_ia():
    """ Plot some Ia SNe from ZTF, using Yuhan's light curves
    and also a IIb """
    ax.scatter(10, -19, marker='D', color='lightblue')
    ax.text(
            10-losx, -19-losy, 'SN Ia', fontsize=12,
            horizontalalignment='right')
    ax.scatter(8, -18.2, marker='D', color='lightblue')
    ax.text(
            8+losx, -18.2+losy, 'SN IIb', fontsize=12)


def plot_ibn():
    """ Ibn SNe: Table 4 in Hosseinzadeh 2017 """
    kwargs = dict(mec='#b2df8a', marker='s', mfc='white', c='#b2df8a')
    akwargs = dict(length_includes_head=True, head_width=0.05, color='#b2df8a')
    dx = -0.5 
    ax.scatter(
            0,0,facecolor='white', 
            edgecolor='#b2df8a', label="Ibn", marker='s') 

    # SN 1999cq
    ax.errorbar(
        3.9, -19.73, yerr=0.10, **kwargs)

    # PTF11rfh
    # Shri says error bar too big
    # ax.errorbar(
    #         7, -20.49, xerr=4.7, yerr=0.99, color=ibn_c, marker=ibn_m, label="Ibn")

    # LSQ12btw
    ax.errorbar(3.8, -19.44, yerr=0.04, **kwargs)
    ax.arrow(3.3, -19.73, dx, 0, **akwargs)

    # PTF12ldy
    ax.errorbar(6.2, -19.20, xerr=2.0, yerr=0.02, **kwargs)

    # iPTF13beo
    ax.errorbar(1.6, -18.57, xerr=0.9, yerr=0.05, **kwargs)

    # LSQ13ccw
    ax.errorbar(4.7, -18.46, xerr=2.1, yerr=0.06, **kwargs)

    # iPTF14aki
    ax.errorbar(7.0, -19.30, xerr=1.9, yerr=0.03, **kwargs)

    # SN 2014bk
    ax.errorbar(9.9, -19.99, yerr=0.65, **kwargs)
    ax.arrow(9.9, -19.99, dx, 0, **akwargs)

    # SN 2015U
    ax.errorbar(8.8, -19.41, xerr=0.9, yerr=0.27, **kwargs)

    # iPTF15ul
    ax.errorbar(3, -20.43, xerr=0.7, yerr=0.26, **kwargs)

    # iPTF15akq
    ax.errorbar(8.3, -18.62, xerr=2.7, yerr=0.31, **kwargs)


def fbot():
    """ The PanSTARRS FBOT (gold and silveR) sample. 
    Taken fro Table 4 of Drout+14 """
    kwargs = dict(mec='#a6cee3', marker='o', mfc='white', c='#a6cee3')
    akwargs = dict(length_includes_head=True, head_width=0.05, color='#a6cee3')
    dx = -0.5

    # for legend
    ax.scatter(
            0,0,facecolor='white', edgecolor='#a6cee3', label="Unclassified") 

    # PS1-10bjp
    ax.errorbar(
            3.4, -18.14, xerr=0.1, yerr=0.11, **kwargs)

    # PS1-11qr
    ax.errorbar(2.9, -19.56, xerr=0.1, yerr=0.08, **kwargs)

    # PS1-11bbq
    ax.errorbar(3.3, -19.73, yerr=0.10, **kwargs)
    ax.arrow(3.3, -19.73, dx, 0, **akwargs)

    # PS1-12bv
    ax.errorbar(2.2, -19.49, yerr=0.07, **kwargs)
    ax.arrow(2.2, -19.49, dx, 0, **akwargs)

    # PS1-12brf
    ax.errorbar(1, -18.43, yerr=0.08, **kwargs)
    ax.arrow(1, -18.43, dx, 0, **akwargs)

    # PS1-13ess
    ax.errorbar(8.9, -18.43, yerr=0.18, **kwargs)
    ax.arrow(8.9, -18.43, dx, 0, **akwargs)

    # Now from DES ... but these are in g-band
    # ax.errorbar(2.6, -19.42, yerr=0.32, c=fbot_c, marker='<')
    # ax.errorbar(6.4, -19.69, yerr=0.15, c=fbot_c, marker='<')
    # ax.errorbar(5.2, -18.98, yerr=0.04, c=fbot_c, marker='<')
    # ax.errorbar(9.2, -19.84, yerr=0.08, c=fbot_c, marker='<')
    # ax.errorbar(6.5, -19.74, yerr=0.06, c=fbot_c, marker='<')
    # ax.errorbar(4.9, -19.76, yerr=0.06, c=fbot_c, marker='<')
    # ax.errorbar(5.8, -18.23, yerr=0.05, c=fbot_c, marker='<')
    # ax.errorbar(7.7, -19.64, yerr=0.02, c=fbot_c, marker='<')
    # ax.errorbar(9.7, -19.46, yerr=0.05, c=fbot_c, marker='<')
    # ax.errorbar(3.8, -18.41, yerr=0.05, c=fbot_c, marker='<')
    # ax.errorbar(4.3, -22.24, yerr=0.07, c=fbot_c, marker='<')
    # ax.errorbar(3.5, -19.84, yerr=0.08, c=fbot_c, marker='<')
    # ax.errorbar(4.5, -20.97, yerr=0.05, c=fbot_c, marker='<')
    # ax.errorbar(5.0, -19.62, yerr=0.05, c=fbot_c, marker='<')

    # Now from Arcavi 2016
    # ax.errorbar(5.87, -20.34, yerr=0.02, c=fbot_c, marker='<')


def kepler():
    """ Plot the Kepler transient """
    ax.scatter(2.2, -18.8, marker='D', c='lightblue'), 
    ax.text(2.2+losx, -18.8+losy, "15K", fontsize=12)



def ptf09uj():
    """ Plot PTF09uj """
    ax.scatter(3.5, -19, marker='D', c='lightblue')
    ax.text(3.5+losx, -19+losy, "09uj", fontsize=12)


def cow():
    """ Plot 18cow """
    x = 1.5
    y = -19.9
    ax.errorbar(x, y, xerr=1, marker='*', c='#1f78b4', markersize=30)
    #ax.arrow(x, y, -0.5, 0, length_includes_head=True, head_width=0.05, 
    #        color='#1f78b4')
    ax.text(x+losx, y+losy, "18cow", fontsize=16)


def asu():
    """ Plot iPTF16asu, this is g-band though """
    x = 4
    y = -20.4
    ax.scatter(x, y, marker='*', c='#1f78b4', s=500)
    #ax.arrow(x, y, -0.5, 0, length_includes_head=True, head_width=0.05, 
    #        color='#1f78b4')
    ax.text(x+losx, y-2*losy, "PTF16asu", fontsize=16, verticalalignment='top',
            horizontalalignment='center')


def koala():
    """ Plot koala"""
    x = 1.5
    y = -20.6
    ax.errorbar(x, y, xerr=0.5, fmt='*', c='#1f78b4', markersize=30)
    #ax.arrow(x, y, -0.5, 0, length_includes_head=True, head_width=0.05, 
    #        color='#1f78b4')
    ax.text(x+losx, y+losy, "ZTF18abvkwla", fontsize=16)


def sn06aj():
    """ Plot SN2006aj """
    x = 6.4/24
    y = -18.4
    ax.scatter(x, y, marker='*', c='#1f78b4', s=800)
    ax.scatter(0, 0, marker='*', c='#1f78b4', s=200, label="Candidate Engines")
    ax.arrow(x, y, -0.5, 0, length_includes_head=True, head_width=0.05, 
            color='#1f78b4')
    ax.text(x+losx, y+losy, "SN2006aj", fontsize=12)


def gep():
    """ Plot gep"""
    x = 3
    y = -19.5
    ax.scatter(x, y, marker='*', c='#1f78b4', s=500, zorder=5)
    #ax.arrow(x, y, -0.5, 0, length_includes_head=True, head_width=0.05, 
    #        color='#1f78b4', zorder=5)
    ax.text(x+losx, y+losy, "ZTF18abukavn", fontsize=16)



plot_ia()
gep()
koala()
cow()
plot_ibn()
asu()
fbot()
kepler()
ptf09uj()
sn06aj()

#plt.axvspan(0 ,5 , facecolor='white', edgecolor='k', lw=3)
#plt.axvline(x=3, ls='--', c='k')

ax.legend(loc='upper right', fontsize=14)
ax.set_xlabel("Rise Time (days)", fontsize=16)
ax.set_ylabel("Peak Magnitude", fontsize=16)
ax.set_xlim(0,10.5)
ax.set_ylim(-21, -18)
ax.invert_yaxis()

plt.tick_params(axis='both', labelsize=14)
plt.tight_layout()

#plt.show()
plt.savefig("M_trise.png", dpi=500)
