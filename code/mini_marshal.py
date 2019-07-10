""" My own mini marshal """

import matplotlib.pyplot as plt
import pandas
import numpy as np
from astropy.io import fits
from penquins import Kowalski
from ztfquery import query, marshal

names = np.loadtxt("/Users/annaho/Dropbox/Projects/Research/Koala/code/keep_field_791_792.txt", dtype=str)

# Initialize the table
print('<!doctype html>')
print('<html>')
print('<head>')
print('<title>SOURCE MARSHAL</title>')
print('</head>')
print('<body>')

print('<font size="5" color="green">Fast Transient Candidates</font>')
print('</br>')
print('<hr>')
print('</br>')

# Connect
# zquery = query.ZTFQuery()
# m = marshal.MarshalAccess()
username = 'ah'
password = 'TetraodonInsists'
s = Kowalski(
        protocol='https', host='kowalski.caltech.edu', port=443,
        verbose=False, username=username, password=password)

# One source per line
hspace = 0

# if there's only one...
if names.size == 1:
    names = [names.tolist()]

for name in names:
    # Print the name
    print('<b><font size="3">%s</font></b>'%(name))
    print('</br>')

    # Light curve was generated during the search
    print('<img src="%s_lc.png" hspace=%s>'%(name,hspace))

    # Diff image cutout was generated during the search
    print('<img src="%s_cutout.png">'%(name))

    # PS1 image cutout was generated during the search
    print('<img src="%s_ps1.png">'%(name))

    # DECaLS image
    print('<img src="%s_ls.jpg">'%(name))

    print('</br>')
    print('</br>')

print('</body>')
print('</html>')
