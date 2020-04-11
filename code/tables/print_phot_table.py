""" Print table of photometry, ground-based and UVOT """

import numpy as np
from math import floor, log10
from astropy.time import Time
from astropy.cosmology import Planck15
from astropy.io import ascii

def round_sig(x, sig=2):
    print(x)
    if x < 0:
        return -round(-x, sig-int(floor(log10(-x)))-1)
    return round(x, sig-int(floor(log10(x)))-1)


def ndec(num):
    dec = str(num).split('.')[-1]
    return len(dec)


d = Planck15.luminosity_distance(z=0.2714).cgs.value

headings = np.array(
        ['Date (MJD)', '$\Delta t$', 'Filter', 'AB Mag'])
label = "opt-phot"
caption = "Optical photometry for ZTF18abvkwla \
from forced photometry on P48 images \citep{Yao2019}. \
Values have not been corrected for Galactic extinction."

# Print the table headers
ncol = len(headings)
colstr = ""
colstr += 'l'
for col in np.arange(ncol-1): colstr+="r"
print(colstr)

colheadstr = ""
for col in np.arange(ncol-1):
    colheadstr += "\colhead{%s} & " %headings[col]
colheadstr += "\colhead{%s}" %headings[-1]

rowstr = ""
for col in np.arange(ncol-1):
    rowstr += "%s & "
rowstr += "%s \\\ \n"

outputf = open("table_%s.txt" %label, "w")
outputf.write("\\begin{deluxetable}{%s} \n" %colstr)
outputf.write("\\tablecaption{%s\label{tab:%s}} \n" %(caption,label))
outputf.write("\\tablewidth{0pt} \n")
outputf.write("\\tablehead{ %s } \n" %colheadstr)
#outputf.write("\\rotate \n")
outputf.write("\\tabletypesize{\scriptsize} \n")
outputf.write("\startdata \n")

data_dir = "/Users/annaho/Dropbox/Projects/Research/Koala/data"
t0 = 58373.407523 # time of the first detection, in MJD

# P48 photometry
f = data_dir + "/ZTF18abvkwla_lct.csv"
dat = ascii.read(f, delimiter=',')
jd = dat['jdobs']
filt = np.array([val[2:3] for val in dat['filter']])
mag = dat['mag']
det = mag > 0
mjd = jd[det]-2400000.5
mag = mag[det]
emag = dat['mag_unc'][det]
ndet = sum(det)
tel = np.array(['P48']*ndet)
filt = filt[det]

# generate dt
dt = mjd-t0

# sort
order = np.argsort(dt)
mjd = mjd[order]
dt = dt[order]
mag = mag[order]
emag = emag[order]
tel = tel[order]
filt = filt[order]

# remove duplicate entries
mjd,ind = np.unique(mjd, return_index=True)
dt = dt[ind]
mag = mag[ind]
emag = emag[ind]
tel = tel[ind]
filt = filt[ind]

# Two preceding upper limits
row = rowstr % (58372.389800000004, -1.02, '$r$', '<21.39')
outputf.write(row)
row = rowstr % (58372.42059999984, -0.99, '$g$', '<21.56')
outputf.write(row)

for ii in np.arange(len(dt)):
    if(dt[ii] < 55):
        mjd_str = '{:<08f}'.format(round_sig(mjd[ii], 11)) # pad with zeros
        dt_str = '{:.2f}'.format(np.round(dt[ii], 2))
        mag_str = '{:.2f}'.format(round_sig(mag[ii], 4))
        emag_str = '{:.2f}'.format(np.round(emag[ii], ndec(mag_str)))
        row = rowstr %(
                mjd_str, dt_str, "$"+filt[ii]+"$", 
                '$%s \pm %s$' %(mag_str, emag_str))
        outputf.write(row)

outputf.write("\enddata \n")
outputf.write("\end{deluxetable} \n")
