""" Check the 1DC LC of the surviving sources (not just programid=3)

Only keep sources that have,
- detections on three consecutive nights
- first detection fainter than peak
- last detection fainter than peak
"""

import numpy as np
import glob
import matplotlib.pyplot as plt
from print_diagnostics import get_lc

search = "recent*field*nostars.txt"
allf = glob.glob(search)

for f in allf:
    field = f.split("_")[1]
    print(field)
    outf = "recent_" + field+"_goodlc.txt"
    if len(glob.glob(outf)) == 0:
        cands = np.loadtxt(f, dtype=str)
        ucands = np.unique(cands)
        if len(ucands) > 0:
            print(len(ucands))
            keep = []
            for name in ucands:
                jd_all,dt_all,mag_all,emag_all,filt_all,pid_all = get_lc(name)

                # one filter
                choose = np.logical_and(pid_all == 3, filt_all == 1)

                # First, check that there are at least three detections that count
                if sum(choose) > 2:
                    jd = jd_all[choose]
                    dt = dt_all[choose]
                    mag = mag_all[choose]
                    emag = emag_all[choose]
                    filt = filt_all[choose]
                    pid = pid_all[choose]
                    # next, calculate the peak mag
                    peak = np.min(mag)

                    # check that the first mag < the peak mag 
                    shape = np.logical_and(
                            mag[0] > peak, mag[-1] > peak)
                    dmag = np.array([mag[0]-peak, mag[-1]-peak])
                    if np.logical_and(shape, max(dmag)>1):
                        keep.append(name)

                # one filter
                choose = np.logical_and(pid_all == 3, filt_all == 2)
                jd = jd_all[choose]
                dt = dt_all[choose]
                mag = mag_all[choose]
                emag = emag_all[choose]
                filt = filt_all[choose]
                pid = pid_all[choose]

                # First, check that there are at least three detections that count
                if sum(choose) > 2:
                    # next, calculate the peak mag
                    peak = np.min(mag)

                    # check that the first mag < the peak mag
                    shape = np.logical_and(
                            mag[0] > peak, mag[-1] > peak)
                    dmag = np.array([mag[0]-peak, mag[-1]-peak])
                    if np.logical_and(shape, max(dmag)>1):
                        keep.append(name)
                 
            keep = np.unique(np.array(keep))
            print(len(keep))
            np.savetxt(outf, keep, fmt='%s')
