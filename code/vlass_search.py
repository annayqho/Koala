import sys
sys.path.append("/Users/annaho/Github/Query_VLASS")
from vlass_search import search_vlass
from ret_radio import get_transients


def run_ret():
    names, ra_raw, dec_raw, dates, z = get_transients()
    ra = []
    dec =[]

    limits = np.zeros(len(names))
    obsdates = np.zeros(len(names))

    for ii,val in enumerate(ra_raw):
        arr = val.split(":")
        ra_temp = arr[0]+"h"+arr[1]+"m"+arr[2]+"s"
        dec_temp = arr[0]+"d"+arr[1]+"m"+arr[2]+"s"

        # Convert dec to degrees and RA to hours
        c = SkyCoord(ra_temp, dec_temp, frame='icrs')
        ra_h = c.ra.hour
        ra.append(ra_h)
        dec_d = c.dec.deg
        dec_append(dec_h)
    search_vlass(names, ra, dec, dates)


def run_maxi():
    inputf = "/Users/annaho/Github/MAXI/MAXI/data/musst_loc.txt"
    dat = Table.read(inputf, format='ascii')
    names = np.array(dat['MUSST_ID'])
    dates = Time(dat['Date'])
    ra = dat['RA']
    dec = dat['Dec']
    c = SkyCoord(ra, dec, frame='icrs', unit='deg')
    ra_h = c.ra.hour
    dec_d = c.dec.deg
    search_vlass(names, ra_h, dec_d, dates)
