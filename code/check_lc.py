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

search = "field*nostars.txt"
allf = glob.glob(search)
f = allf[0]
field = f.split("_")[0]
outf = field+"_goodlc.txt"
cands = np.loadtxt(f, dtype=str)
print(len(cands))

keep = []
for name in cands:
    jd,dt,mag,emag,filt,pid = get_lc(name)
    choose = pid == 3
    jd = jd[choose]
    dt = dt[choose]
    mag =mag[choose]
    emag = emag[choose]
    filt = filt[choose]

    # First, check that there are at least three detections that count
    ndet = sum(pid==3)

    if ndet > 2:
        # next, calculate the peak mag
        peak = np.min(mag)

        # check that the first mag < the peak mag
        shape = np.logical_and(mag[0] > peak, mag[-1] > peak)
        if shape:
            keep.append(name)

np.savetxt(outf, keep, fmt='%s')
