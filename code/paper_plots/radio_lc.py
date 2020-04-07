import numpy as np
import matplotlib
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MultipleLocator,NullFormatter
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)


def plot_band(ax, t, f, ef, c, marker, label, limit=None):
    """
    Plot a single band of radio data

    Parameters
    ----------
    t: time in observer frame
    f: flux in observerframe, mJy
    ef: uncertainty on flux in observer frame, mJy
    c: color of points
    marker: marker of points
    label: label for legend
    limit: optional, an array indicating which points are limits
    """
    ax.errorbar(
            t, f, ef, mec='k', mfc=c, c='k', marker=marker, label=label, ms=8)
    if limit is not None:
        xval = t[limit][0]
        yval = f[limit][0]
        ax.arrow(
                xval, yval, 0, -yval/5, color='k', length_includes_head=True, 
                head_width=xval/13, head_length=yval/10)
    if len(t)>1:
        ax.plot(t, f, c=c, marker=marker)


def plot_sed(ax):
    axins = inset_axes(
            ax, 2, 1, loc=3,
            bbox_to_anchor=(0.15,0.15),
            bbox_transform=ax.transAxes)
    freq = np.array([10,6,3,1.5,0.65])
    flux = np.array([0.045, 0.089, 0.067, 0.135, 0.105])
    eflux = np.array([0.003,0.003,0.005,0.007,0])
    marker = np.array(['s', 'D', 'o', '^', 'o'])
    c = np.array(['k', '#f98e09', 'white', '#57106e', 'k'])
    for ii,m in enumerate(marker):
        axins.errorbar(freq[ii], flux[ii], eflux[ii], 
                fmt=m, mfc=c[ii], mec='k', c='k')
    axins.arrow(freq[-1], flux[-1], 0, -0.01, 
            color='k', head_width=0.7, head_length=0.01)
    axins.tick_params(axis='both', labelsize=12)
    axins.set_xlabel("Observed Frequency (GHz)", fontsize=12)
    axins.text(0.9,0.9,r'$\Delta t\approx350$\,d', transform=axins.transAxes,
            fontsize=11, horizontalalignment='right', verticalalignment='top')


if __name__=="__main__":
    fig,ax = plt.subplots(1,1,figsize=(6,4))

    # X-band
    plot_band(ax, np.array([81, 310, 352, 396]), 
            np.array([0.364, 0.061, 0.045, 0.031]), 
            np.array([0.006, 0.003, 0.003, 0.003]), 'k', 's', "VLA: 10 GHz")

    # C-band
    plot_band(ax, np.array([343,397]), np.array([0.089,0.033]), 
            np.array([0.003,0.003]), '#f98e09', 'D', label="VLA: 6 GHz")

    # S-band
    plot_band(ax,
            np.array([188,346]), np.array([0.134,0.067]), np.array([0,0.005]), 
            'white', 'o', "VLA: 3 GHz", limit=np.array([True,False]))

    # L-band
    plot_band(ax,
            np.array([351]), np.array([0.135]), np.array([0.007]),
            '#57106e', '^', "VLA: 1.5 GHz")

    # GMRT point
    plot_band(ax,
            np.array([364]), np.array([0.105]), np.array([0]),
            'k', 'o', "GMRT: 650 MHz", limit=np.array([True]))

    plot_sed(ax)

    ax.set_ylabel(r"$f_\nu$ (mJy)", fontsize=16)
    ax.set_xlabel("Observer-Frame Days Since 2018 Sep 11.42", fontsize=16)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_yticks([0.1,0.3,0.05])
    ax.set_xticks([100,200,300])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.tick_params(axis='both', labelsize=16)
    ax.legend(fontsize=10, loc='upper right', ncol=2)
    ax.xaxis.set_minor_formatter(NullFormatter())

    fig.tight_layout()
    plt.savefig("radio_lc.png", dpi=300)
    #plt.show()
