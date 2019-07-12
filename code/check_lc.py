""" Check the 1DC LC of the surviving sources (not just programid=3)

Only keep sources that have,
- detections on three consecutive nights
- first detection fainter than peak
- last detection fainter than peak
"""

import numpy as np
import glob
from print_diagnostics import get_lc

search = "field*nostars.txt"
allf = glob.glob(search)
f = allf[0]
cands = np.loadtxt(f, dtype=str)
name = cands[0]

jd,dt,mag,emag,filt = get_lc(name)
