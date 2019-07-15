import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from ztfquery import query,marshal
from print_diagnostics import get_lc


def rise(jd_all, dt_all, mag_all, emag_all, pid_all, filt_all, filt=1):
    """ Measure the fade time. Default filter is g. """
    choose = np.logical_and(pid_all != 2, filt_all == filt)
    jd = jd_all[choose]
    dt = dt_all[choose]
    mag = mag_all[choose]
    emag = emag_all[choose]
    pid = pid_all[choose]

    if sum(choose)>0:
        ind = np.argmin(mag)
        mpeak = mag[ind]
        tpeak = dt[ind] 

        dt_fit = dt[0:ind+1]
        mag_fit = mag[0:ind+1]

        if len(dt_fit) > 1:
            dmag = mag_fit[0]-mag_fit[-1]
            if np.logical_and(dmag < 1, dt_fit[-1]-dt_fit[0]>=7):
                # if it hasn't reached dmag=1 
                # but you can still tell that it's bad
                return 2
            else:
                order = np.argsort(mag_fit)
                dt_fit = dt_fit[order]
                mag_fit = mag_fit[order]
                tval = np.interp(mpeak+1, mag_fit, dt_fit, left=-999, right=-999)
                if tval != -999:
                    trise = tpeak - tval
                    if trise < 7:
                        return 1
                    else:
                        return 2
    return 0


def fade(jd_all, dt_all, mag_all, emag_all, pid_all, filt_all, filt=1):
    """ Measure the fade time. Default filter is g. """
    choose = np.logical_and(pid_all != 2, filt_all == 1)
    if sum(choose)>0:
        jd = jd_all[choose]
        dt = dt_all[choose]
        mag = mag_all[choose]
        emag = emag_all[choose]

        ind = np.argmin(mag)
        mpeak = mag[ind]
        tpeak = dt[ind] 

        dt_fit = dt[ind:]
        mag_fit = mag[ind:]

        if len(dt_fit) > 1:
            dmag = mag_fit[-1]-mag_fit[0]
            if np.logical_and(dmag < 1, dt_fit[-1]-dt_fit[0]>=8):
                return 2
            else:
                order = np.argsort(mag_fit)
                dt_fit = dt_fit[order]
                mag_fit = mag_fit[order]

                tval = np.interp(mpeak+1, mag[ind:][order], dt[ind:][order], left=-999, right=-999)
                if tval != -999:
                    tfade = tval - tpeak
                    print(tfade)
                    if tfade < 8:
                        return 1
                    else:
                        return 2
    return 0


def load_data():
    m = marshal.MarshalAccess()    
    
    cands = np.loadtxt("passed_lc_check.txt", dtype=str)
    m = marshal.MarshalAccess()
    tab = m.get_target_data(cands)
    
    names = tab['name']
    classes = tab['classification']
    classes = np.array(classes)
    np.savetxt("passed_lc_check.txt", np.array([names,classes]).T, fmt='%s',
            delimiter=',')


if __name__=="__main__":
    names = np.loadtxt("real.txt", dtype=str)
    keep = []

    names = ['ZTF19aapfmki']
    for ii,sn in enumerate(names):
        # 0: no info, 1: short, 2: long
        print(ii)
        print(sn)
        jd_all,dt_all,mag_all,emag_all,filt_all,pid_all = get_lc(sn)
        trise_g = rise(jd_all, dt_all, mag_all, emag_all, pid_all, filt_all, 1)
        trise_r = rise(jd_all, dt_all, mag_all, emag_all, pid_all, filt_all, 2)
        tfade_g = rise(jd_all, dt_all, mag_all, emag_all, pid_all, filt_all, 1)
        tfade_r = fade(jd_all, dt_all, mag_all, emag_all, pid_all, filt_all, 2)

        # You would only keep it if there is at least one 1
        if np.logical_or.reduce((trise_g==1, trise_r==1, tfade_g==1, tfade_r==1)):
            # and no 2's
            if np.logical_and.reduce((trise_g!=2, trise_r!=2, tfade_g!=2, tfade_r!=2)):
                print("Keeping this one")
                keep.append(sn)

    np.savetxt("real_lc_check.txt", keep, fmt='%s')
