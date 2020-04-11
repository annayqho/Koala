import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.io import ascii
from astropy.cosmology import Planck15
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

z = 0.2714
rcol = '#d95f02'
gcol = '#1b9e77'
t0 = 58373.407523+2400000.5

# Full light curve
def full(ax):
    dat = ascii.read("../../data/koala_lc.txt")
    dt = np.array(dat['col1'])-t0
    flux = np.array(dat['col2'])
    eflux = np.array(dat['col3'])
    filt = np.array(dat['col4'])

    choose = filt == 'r'
    ax.errorbar(dt[choose]/(1+z),flux[choose],eflux[choose],mec=rcol,fmt='o',
            mfc='white',c=rcol,
            label=r"ZTF $r$ ($\lambda_\mathrm{rest}=4950$\AA)")

    choose = filt == 'g'
    ax.errorbar(dt[choose]/(1+z),flux[choose],eflux[choose],c=gcol,fmt='s',
            label=r"ZTF $g$ ($\lambda_\mathrm{rest}=3820$\AA)")

    # Plot the radioactively powered LC
    mni_vals = [0.1, 0.2, 0.3]
    for mni in mni_vals:
        tmod = np.linspace(5,40,10000)
        lum = 2E43 * (mni) * \
            (3.9*np.exp(-tmod/8.8)+0.678*(np.exp(-tmod/113.6)-np.exp(-tmod/8.8)))
        lam = 4950*1E-8 # cm
        nu = 3E10/lam # Hz
        dcm = Planck15.luminosity_distance(z=0.2714).cgs.value
        flam = lum/nu/(4*np.pi*dcm**2)
        fmcmc = flam*10**(29.44)
        ax.plot(tmod, fmcmc, c='k', ls='--', lw=0.5)
        if mni==0.3:
            ax.text(
                    5, fmcmc[0], "$M_\mathrm{Ni}=$"+str(mni)+r"$M_\odot$", 
                    fontsize=12, 
                    horizontalalignment='right', verticalalignment='center')
        else:
            ax.text(
                    5, fmcmc[0], str(mni)+r"$M_\odot$", fontsize=12, 
                    horizontalalignment='right', verticalalignment='center')
    ax.axhline(y=0, c='lightgrey')

    ax.set_xlabel("Rest-frame days since first detection", fontsize=16)
    ax.set_ylabel("Flux (ZP=25)", fontsize=16)
    ax.tick_params(axis='both', labelsize=16)
    ax.set_ylim(-50,170)

    ax.legend(fontsize=13, loc='lower left', ncol=1)

     

# Magnitude space
def mag(ax):
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
            dt, gmag, yerr=egmag, fmt='s', c=gcol, lw=0.5, ms=8, 
            label=r"ZTF $g$ ($\lambda_\mathrm{rest}=3820$\AA)")

    # Plot the g-band upper limits
    choose = np.logical_and(filts=='g', dat['mag'] == 0)
    jd = dat['jdobs'][choose]
    dt = (jd-t0)/(1+z)
    lims = dat['limmag'][choose] - 0.167 # extinction in SDSS g
    ax.scatter(
            dt, lims, edgecolor=gcol, facecolor='white', marker='s', label=None)
    for ii,dt_val in enumerate(dt):
        ax.arrow(dt_val, lims[ii], 0, +0.1, color=gcol, 
                length_includes_head=True,
                head_width=0.2, head_length=0.05, label=None)
     
    # Plot the r-band light curve
    choose = np.logical_and(filts=='r', dat['mag'] > 0)
    jd = dat['jdobs'][choose]
    dt = (jd-t0)/(1+z)
    rmag = dat['mag'][choose] - 0.115
    ermag = dat['mag_unc'][choose]
    ax.errorbar(
            dt, rmag, yerr=ermag, fmt='o', mec=rcol, mfc=rcol, c=rcol, 
            lw=0.5, ms=8,
            label=r"ZTF $r$ ($\lambda_\mathrm{rest}=4950$\AA)")

    # Plot the r-band upper limits
    choose = np.logical_and(filts=='r', dat['mag'] == 0)
    jd = dat['jdobs'][choose]
    dt = (jd-t0)/(1+z)
    lims = dat['limmag'][choose] - 0.115
    ax.scatter(dt, lims, edgecolor=rcol, facecolor='white', 
            marker='o', label=None)
    for ii,dt_val in enumerate(dt):
        ax.arrow(dt_val, lims[ii], 0, +0.1, edgecolor=rcol, lw=0.5,
                facecolor=rcol, length_includes_head=True,
                head_width=0.2, head_length=0.05, label=None)

    # Plot the 18cow light curves
    # scale to the distance of Koala
    dmag = Planck15.distmod(z=0.2714).value-Planck15.distmod(0.0141).value

    # what filter is closest to 4950 rest-frame?
    # for LT, it's in g-band
    t_g = np.array([58287.1500, 58288.2020, 58289.0234, 58289.1889, 58289.9147,
        58290.1082, 58290.2189, 58290.2922, 58290.9029, 58291.0145, 58291.1742])
    t0_cow = t_g[0] - 2
    dt_g = t_g - t0_cow
    mag_g = np.array([13.40, 13.65, 14.10, 14.18, 14.48, 
        14.57, 14.63, 14.70, 14.94, 14.99, 15.09]) - 0.287 # extinction in g
    ax.plot(dt_g, mag_g + dmag, c=rcol, lw=0.5)
    ax.text(4.3, 21.0, "18cow $g$", fontsize=11, rotation=-55)

    # Connect the LC to the last upper limit 
    # doesn't look good to do this...
    #ax.plot(
    #        [58284.1300-t0_cow, dt_g[0]], 
    #        [18.90 + dmag, mag_g[0]+dmag], c='k', lw=0.5, ls='--')
    # instead, just indicate the last non-detection
    ax.plot([-1, -1], [21.5, 21.8], c=rcol, lw=2)
    ax.text(-0.9, 21.68, 'time of last 18cow upper limit', fontsize=11)

    # what filter is closest to 3820 rest-frame? I think u-band
    # but it's sort of between u-band and b-band
    # for LT, it's pretty squarely in u-band
    t_u = np.array([58289.0227, 58289.1963, 58289.6434, 58289.9154, 58290.2243,
        58290.2267, 58290.6697, 58290.9036, 58291.1809])
    mag_u = np.array(
            [13.97, 14.03, 14.24, 14.31, 14.43, 14.46, 14.63, 14.66, 14.83])
    ax.plot(t_u-t0_cow, mag_u - 0.368 + dmag, c=gcol, lw=0.5, ls='-')
    ax.text(5, 21.05, "18cow $u$", fontsize=11, rotation=-55)

    # ni-powered light curve
    mni_vals = [0.3, 0.4]
    tmod = np.linspace(3,8,100)
    for mni in mni_vals:
        tmod = np.linspace(5,20,10000)
        lum = 2E43 * (mni) * \
            (3.9*np.exp(-tmod/8.8)+0.678*(np.exp(-tmod/113.6)-np.exp(-tmod/8.8)))
        lam = 4950*1E-8 # cm
        nu = 3E10/lam # Hz
        dcm = Planck15.luminosity_distance(z=0.2714).cgs.value
        flux= lum/nu/(4*np.pi*dcm**2)
        m = -2.5*np.log10(flux)-48.6
        ax.plot(tmod,m,c='k',ls='--',lw=0.5)
        ax.text(
                8.4, m[0]*1.009, str(mni)+r"$M_\odot$", fontsize=12, 
                horizontalalignment='left', verticalalignment='top')

    # Format the panel
    ax.set_ylim(21.7, 19.2)
    ax2 = ax.twinx()
    ax2.set_ylabel(
            r"Absolute Mag",
            fontsize=14, rotation=270, labelpad=15.0)
    y_f = lambda y_i: y_i-Planck15.distmod(z=0.2714).value
    ymin, ymax = ax.get_ylim()
    print(ymin)
    print(y_f(ymin))
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.plot([],[])
    ax2.tick_params(axis='both', labelsize=14)

    ax.legend(fontsize=13, loc='upper right', ncol=1)

    ax.text(0.9, 19.32, 'S', fontsize=12)

    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlim(-2.2,20)
    ax.set_xlabel("Rest-frame days since first detection", fontsize=16)
    ax.set_ylabel("Apparent Mag", fontsize=16)
    ax.tick_params(axis='both', labelsize=16)


if __name__=="__main__":
    fig,ax = plt.subplots(1,1, figsize=(11,5))
    mag(ax)
    plt.tight_layout()
    #plt.show()
    plt.savefig("lc.png", dpi=300)
