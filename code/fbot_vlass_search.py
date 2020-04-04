import numpy as np
from astropy.coordinates import SkyCoord
import sys
sys.path.append("/Users/annaho/Github/Query_VLASS")
from vlass_search import run_search
from ret_radio import get_transients


def run_ret():
    names, ra_raw, dec_raw, dates, z = get_transients()
    ra = []
    dec =[]

    limits = np.zeros(len(names))
    dts = np.zeros(len(names))

    for ii,val in enumerate(ra_raw):
        arr = val.split(":")
        ra_temp = arr[0]+"h"+arr[1]+"m"+arr[2]+"s"
        arr2 = dec_raw[ii].split(":")
        dec_temp = arr2[0]+"d"+arr2[1]+"m"+arr2[2]+"s"

        # Convert dec to degrees and RA to hours
        c = SkyCoord(ra_temp, dec_temp, frame='icrs')
        lim,obsdate = search_vlass(names[ii], c, dates[ii])
        limits[ii] = lim
        dts[ii] = obsdate-dates[ii]
    return names, z, ra_raw, dec_raw, dts, limits


if __name__=="__main__":
    names, z, ra_raw, dec_raw, dts, limits = run_ret()
