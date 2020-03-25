""" Measure the effective volume of the survey by taking the seeing
from night to night """

import numpy as np
from ztfquery import query
from astropy.io import ascii
from astropy.time import Time
from astropy.cosmology import Planck15


def get_seeing(zquery):
    dat = ascii.read("../data/fieldnights_1DC.txt", "w")
    f = dat['col1'] # fields 
    fn = dat['col2'] # column-nights

    outputf = open("seeing.txt", "w")

    for ii,field in enumerate(f):
        startjd = Time(fn[ii]).jd
        endjd = startjd+1
        query = "field=%s AND obsjd BETWEEN %s AND %s" %(field,startjd,endjd)
        zquery.load_metadata(sql_query=query)
        tab = zquery.metatable
        maglims = tab['maglimit'].values
        if len(maglims) == 0:
            print(field,fn[ii])
        else:
            medlim = np.median(tab['maglimit'].values)
            print(str(ii) + "," + str(medlim))
            outputf.write(str(medlim) + '\n')
    outputf.close()


if __name__=="__main__":
    #zquery = query.ZTFQuery()
    #get_seeing(zquery)

    # generate a grid of limiting magnitude and corresponding redshift
    maglim_grid = np.arange(16, 24, 1)
    z = np.array([0.035, 0.0546, 0.0847, 0.1302, 0.198, 0.298, 0.441, 0.647])
    # this is linear in the log(z) vs. limiting magnitude space

    maglim = np.loadtxt("seeing.txt", delimiter=',')[:,1]

    # for each limiting magnitude, you have an effective volume for a transient
    # with a light curve like AT2018cow. Let's say that you would need to get
    # to an absolute magnitude of -20.
    zlim = np.zeros(len(maglim))
    dlim = np.zeros(len(maglim))
    for ii,lim in enumerate(maglim):
        zval = 10**(np.interp(lim, maglim_grid, np.log10(z)))
        zlim[ii] = zval
        dlim[ii] = Planck15.luminosity_distance(zval).value
