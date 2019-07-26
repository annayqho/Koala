""" Plot absolute magnitude vs. rise time

I'm interested in knowing what peaks between -18 and -21 mag
with rise times between 1 day and 10 days
"""


import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
from astropy.cosmology import Planck15


fig, ax = plt.subplots(1,1, figsize=(6,4))


def plot_ibn():
    """ Ibn SNe: Table 4 in Hosseinzadeh 2017 """
    ibn_c = 'lightgrey'
    ibn_m = 's'

    # SN 1999cq
    ax.errorbar(3.9, -19.73, yerr=0.10, color=ibn_c, marker='<')

    # PTF11rfh
    ax.errorbar(
            7, -20.49, xerr=4.7, yerr=0.99, color=ibn_c, marker=ibn_m, label="Ibn")

    # LSQ12btw
    ax.errorbar(3.8, -19.44, yerr=0.04, color=ibn_c, marker='<')

    # PTF12ldy
    ax.errorbar(6.2, -19.20, xerr=2.0, yerr=0.02, color=ibn_c, marker=ibn_m)

    # iPTF13beo
    ax.errorbar(1.6, -18.57, xerr=0.9, yerr=0.05, color=ibn_c, marker=ibn_m)

    # LSQ13ccw
    ax.errorbar(4.7, -18.46, xerr=2.1, yerr=0.06, color=ibn_c, marker=ibn_m)

    # iPTF14aki
    ax.errorbar(7.0, -19.30, xerr=1.9, yerr=0.03, color=ibn_c, marker=ibn_m)

    # SN 2014bk
    ax.errorbar(9.9, -19.99, yerr=0.65, color=ibn_c, marker='<')

    # SN 2015U
    ax.errorbar(8.8, -19.41, xerr=0.9, yerr=0.27, color=ibn_c, marker=ibn_m)

    # iPTF15ul
    ax.errorbar(3, -20.43, xerr=0.7, yerr=0.26, color=ibn_c, marker=ibn_m)

    # iPTF15akq
    ax.errorbar(8.3, -18.62, xerr=2.7, yerr=0.31, color=ibn_c, marker=ibn_m)


def fbot():
    """ The PanSTARRS FBOT (gold and silveR) sample. 
    Taken fro Table 4 of Drout+14 """
    fbot_c = 'orange' 
    fbot_m = 'o'

    # PS1-10bjp
    ax.errorbar(
            3.4, -18.14, xerr=0.1, yerr=0.11, 
            c=fbot_c, marker=fbot_m, label="PS1")

    # PS1-11qr
    ax.errorbar(2.9, -19.56, xerr=0.1, yerr=0.08, c=fbot_c, marker=fbot_m)

    # PS1-11bbq
    ax.errorbar(3.3, -19.73, yerr=0.10, c=fbot_c, marker='<')

    # PS1-12bv
    ax.errorbar(2.2, -19.49, yerr=0.07, c=fbot_c, marker='<')

    # PS1-12brf
    ax.errorbar(1, -18.43, yerr=0.08, c=fbot_c, marker='<')

    # PS1-13ess
    ax.errorbar(8.9, -18.43, yerr=0.18, c=fbot_c, marker='<')


def kepler():
    """ Plot the Kepler transient """
    ax.scatter(2.2, -18.8, marker='D', c='lightblue'), 
    ax.text(2.2, -18.8, "15K", fontsize=12)



def ptf09uj():
    """ Plot PTF09uj """
    ax.scatter(3.5, -19, marker='D', c='lightblue')
    ax.text(3.5, -19, "09uj", fontsize=12)


def cow():
    """ Plot 18cow """
    ax.scatter(1.5, -19.9, marker='D', c='lightblue')
    ax.text(1.5, -19.9, "18cow", fontsize=12)


def koala():
    """ Plot koala"""
    ax.scatter(2, -20.6, marker='<', c='black', s=50)
    ax.text(2, -20.6, "Koala (ZTF)", fontsize=14)


def gep():
    """ Plot gep"""
    ax.scatter(3, -19.5, marker='<', c='black', s=50)
    ax.text(3, -19.5, "18gep (ZTF)", fontsize=14)


gep()
koala()
cow()
plot_ibn()
fbot()
kepler()
ptf09uj()
ax.legend(loc='upper right', fontsize=14)
ax.set_xlabel("Rise Time (days)", fontsize=16)
ax.set_ylabel("Peak Magnitude (r-band)", fontsize=16)
ax.set_xlim(0,10)
ax.set_ylim(-21, -18)
ax.invert_yaxis()

plt.tick_params(axis='both', labelsize=14)
plt.tight_layout()

plt.savefig("M_trise.png", dpi=500)
