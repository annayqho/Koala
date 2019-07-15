""" Plot the rise timescale and the fade timescale

Say, the time to rise by 1 mag to peak, and the time to fade
by 1 mag from peak
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from ztfquery import query,marshal

m = marshal.MarshalAccess()
jdstart = 2458366
jdend = 2458453


def get_peak(lc_dict, use_filt):
    jds = lc_dict['jdobs'].values
    mags = lc_dict['mag'].values
    isdiffpos = lc_dict['isdiffpos'].values
    filt = lc_dict['filter'].values

    choose = np.logical_and(isdiffpos, filt==use_filt)
    index = np.argmin(mags[choose])

    return jds[choose][np.argmin(mags[choose])], np.min(mags[choose])


def get_rise(lc_dict):
    # interpolate detections in a single band.
    jds = lc_dict['jdobs'].values
    mags = lc_dict['mag'].values
    filt = lc_dict['filter'].values
    isdiffpos = lc_dict['isdiffpos'].values

    rises = []

    jd_peak, m_peak = get_peak(lc_dict, 'g')
    mbase = m_peak+1
    choose = np.logical_and.reduce(
            (mags<99, jds>jdstart, jds<jd_peak, isdiffpos, filt=='g'))

    if sum(choose) > 1:
        # have to sort for np.interp to work
        order = np.argsort(mags[choose])
        startt = np.interp(
                mbase, mags[choose][order], jds[choose][order], 
                left=-99, right=-99)
        if startt > 0:
            rise = jd_peak - startt
            rises.append(rise)

    jd_peak, m_peak = get_peak(lc_dict, 'r')
    mbase = m_peak+1
    choose = np.logical_and.reduce(
            (mags<99, jds>jdstart, jds<jd_peak, isdiffpos, filt=='r'))
    if sum(choose) > 1:
        # have to sort for np.interp to work
        order = np.argsort(mags[choose])
        startt = np.interp(
                mbase, mags[choose][order], jds[choose][order], 
                left=-99, right=-99)
        if startt > 0:
            rise = jd_peak - startt
            rises.append(rise)
     
    if len(rises) == 0:
        return -99
    elif len(rises) == 1:
        return rises[0]
    else:
        return np.min(np.array(rises))


def get_fade(lc_dict):
    # interpolate detections in a single band.
    jds = lc_dict['jdobs'].values
    mags = lc_dict['mag'].values
    filt = lc_dict['filter'].values
    isdiffpos = lc_dict['isdiffpos'].values

    fades = []

    jd_peak, m_peak = get_peak(lc_dict, 'g')
    mbase = m_peak+1
    choose = np.logical_and.reduce(
            (mags<99, jds>jd_peak, jds<jdend, isdiffpos, filt=='g'))

    if sum(choose) > 1:
        endt = np.interp(
                mbase, mags[choose], jds[choose], 
                left=-99, right=-99)
        if endt > 0:
            fade = endt - jd_peak
            fades.append(fade)
    
    jd_peak, m_peak = get_peak(lc_dict, 'r')
    mbase = m_peak+1
    choose = np.logical_and.reduce(
            (mags<99, jds>jd_peak, jds<jdend, isdiffpos, filt=='r'))

    if sum(choose) > 1:
        # have to sort for np.interp to work
        endt = np.interp(
                mbase, mags[choose], jds[choose], 
                left=-99, right=-99)
        if endt > 0:
            fade = endt - jd_peak
            fades.append(fade)

    if len(fades) == 0:
        return -99
    elif len(fades) == 1:
        return fades[0]
    else:
        return np.min(np.array(fades))


if __name__=="__main__":
    rises_all = []
    fades_all = []

    fig,axarr = plt.subplots(1,2, figsize=(8,4))

    cands = np.loadtxt("toscan_byeye.txt", dtype=str)
    m = marshal.MarshalAccess()
    m.get_target_classification(cands)
    classes = []
    for cand in cands:
        classes.append(m.get_target_classification(cand).values[0])
    classes = np.array(classes)
    
    # Choose the Ia SNe and calculate their trise and tfade
    isia = np.array(['Ia' in val for val in classes])
    snia = cands[isia]

    for sn in snia:
        print(cand)
        # get the LC from the marshal
        marshal.download_lightcurve(cand)
        lc_dict = marshal.get_local_lightcurves(cand)
        rise = get_rise(lc_dict)
        fade = get_fade(lc_dict)
        print(rise)
        print(fade)
        rises_all.append(rise)
        fades_all.append(fade)

    ax = axarr[0]
    rises_all = np.array(rises_all)
    ax.hist(rises_all[rises_all>0], histtype='step', color='k')
    ax.axvline(x=2, c='k', lw=2, label="FBOT")
    ax.axvline(x=9, c='r', lw=2, label="Ia SN")
    ax.axvline(x=11, c='r', lw=2)
    ax.tick_params(labelsize=14, axis='both')

    ax.set_xlabel("Rise Time in Days", fontsize=16)
    ax.set_ylabel("Number of Objects", fontsize=16)

    # ax = axarr[1]
    # fades_all = np.array(fades_all)
    # ax.hist(fades_all[fades_all>0], histtype='step', color='k')

    # ax.axvline(x=4.5, c='k', lw=2)
    # ax.axvline(x=18, c='r', lw=2, label="Ia SN")
    # ax.axvline(x=2, c='k', lw=2, label="FBOT")

    # ax.set_xlabel("Fade Time in Days", fontsize=16)
    # ax.set_ylabel("Number of Objects", fontsize=16)
    # ax.tick_params(labelsize=14, axis='both')

    # plt.tight_layout()
    # plt.show()
