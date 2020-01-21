from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)

def plot_not():
    fig,ax = plt.subplots(1,1,figsize=(10,10))

    dat = np.loadtxt("15ul_20150310_FTN_v7.ascii")
    plt.step(dat[:,0][1000:], dat[:,1][1000:]/1E-16+50, where='mid', c='k', lw=0.5)
    plt.text(10000, 56, "FTN, -2", fontsize=12, horizontalalignment='right')

    dat = np.loadtxt("15ul_20150311_NOT_v3.ascii")
    plt.step(dat[:,0], (dat[:,1]/1E-16)+45, where='mid', c='k', lw=0.5)
    plt.text(9000, 50, "NOT, -1", fontsize=12, horizontalalignment='right')

    dat = np.loadtxt("15ul_20150312_NOT_v1.ascii")
    plt.step(dat[:,0], (dat[:,1]/1E-16)+40, where='mid', c='k', lw=0.5)
    plt.text(9000, 45, "NOT, 0", fontsize=12, horizontalalignment='right')

    #dat = np.loadtxt("15ul_20150321_Gemini_N_v2.ascii")
    #plt.step(dat[:,0], (dat[:,1]/1E-16)+40, where='mid', c='k', lw=0.5)
     
    dat = np.loadtxt("15ul_20150328_NOT_v1.ascii")
    plt.step(dat[:,0], (dat[:,1]/1E-16)+35, where='mid', c='k', lw=0.5)
    plt.text(9000, 41, "NOT, 16", fontsize=12, horizontalalignment='right')

    plt.ylim(35,63)
    plt.gca().get_yaxis().set_ticks([])
    plt.tick_params(axis='both', labelsize=14)
    plt.yscale('log')
    plt.xlabel("Wavelength (AA)", fontsize=16)
    plt.ylabel("log Flux + offset", fontsize=16)
    plt.tight_layout()
    plt.savefig("NOT_spec_sequence.png", dpi=200)


def plot_gmos():
    """ plot the GMOS and host galaxy spectra """
    fig,ax = plt.subplots(1,1,figsize=(10,5))

    dat = np.loadtxt("15ul_20150321_Gemini_N_v2.ascii")
    plt.step(
            dat[:,0], dat[:,1]/1E-18, 
            where='mid', c='k', lw=0.5)
    plt.text(8500, 125, "GMOS-N, +9", fontsize=12, horizontalalignment='right')

    dat = fits.open("host_sdss_spec.fits")
    f = dat[1].data['flux']
    wl = 10**(dat[1].data['loglam'])
    plt.step(
            wl, f+70,
            where='mid', c='lightblue', lw=0.5)
    #plt.text(9000, 51.2, "GMOS-N, +9", fontsize=12, horizontalalignment='right')


    plt.ylim(80,240)
    plt.xlim(3800,9000)
    plt.gca().get_yaxis().set_ticks([])
    plt.tick_params(axis='both', labelsize=14)
    plt.yscale('log')
    plt.xlabel("Wavelength (AA)", fontsize=16)
    plt.ylabel("Flux + offset", fontsize=16)
    plt.tight_layout()
    #plt.savefig("GMOS_spec.png", dpi=200)
    plt.show()


def plot_halpha():
    fig,ax = plt.subplots(2,2,figsize=(8,5))

    dat = np.loadtxt("15ul_20150311_NOT_v3.ascii")
    ax[0,0].step(dat[:,0], (dat[:,1]/1E-16)+45, where='mid', c='k', lw=0.5)
    ax[0,0].text(9000, 50, "NOT, -1", fontsize=12, horizontalalignment='right')
    ax[0,0].set_ylim(48,55)
    ax[0,0].text(0.99,0.9,'NOT -1, $R\sim700$',transform=ax[0,0].transAxes,
            horizontalalignment='right')

    dat = np.loadtxt("15ul_20150312_NOT_v1.ascii")
    ax[0,1].step(dat[:,0], (dat[:,1]/1E-16)+40, where='mid', c='k', lw=0.5)
    ax[0,1].text(9000, 45, "NOT, 0", fontsize=12, horizontalalignment='right')
    ax[0,1].set_ylim(43, 52)
    ax[0,1].text(0.99,0.9,'NOT 0, $R\sim700$',transform=ax[0,1].transAxes,
            horizontalalignment='right')

    dat = np.loadtxt("15ul_20150321_Gemini_N_v2.ascii")
    ax[1,0].step(dat[:,0], (dat[:,1]/1E-16)+40, where='mid', c='k', lw=0.5)
    ax[1,0].set_ylim(41, 42.5)
    ax[1,0].set_xlabel("Wavelength (AA)", fontsize=12)
    ax[1,0].text(0.99,0.9,'GMOS +9, $R\sim2000$',transform=ax[1,0].transAxes,
            horizontalalignment='right')
     
    dat = fits.open("host_sdss_spec.fits")
    f = dat[1].data['flux']
    wl = 10**(dat[1].data['loglam'])
    ax[1,1].step(
            wl, f+70, where='mid', c='black', lw=0.5)
    ax[1,1].set_xlabel("Wavelength (AA)", fontsize=12)
    ax[1,1].text(0.99,0.9,'Host eBOSS $R\sim2000$',transform=ax[1,1].transAxes,
            horizontalalignment='right')
    
    for panel in ax.reshape(-1):
        panel.set_xlim(6900,7100)
        panel.set_yticks([])
        panel.tick_params(axis='x', labelsize=12)

    plt.tight_layout()
    plt.show()
     
if __name__=="__main__":
    plot_halpha()
