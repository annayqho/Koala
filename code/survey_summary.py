""" Use ZTFQuery to get basic info about the nightly cadence survey """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.time import Time
from sqlalchemy import create_engine
import pandas as pd

# These are all the 1DC observations
path = '/Users/annaho/Dropbox/Projects/Research/Koala/docs/anna.db'
engine = create_engine(f'sqlite:///{path}')

df = pd.read_sql('Summary', engine)

# Eric says: note that the limiting mags are not meaningful

# Load in the relevant info
mjd = df['expMJD']
fields = df['fieldID']
filt = df['filter']

# Use g-band and Field ID 554 (field with The Koala)
choose = np.logical_and(fields == 554, filt=='g')
g_obsdates = mjd[choose].values
g_dt = g_obsdates[1:]-g_obsdates[0:-1]

choose = np.logical_and(fields == 554, filt=='r')
r_obsdates = mjd[choose].values
r_dt = r_obsdates[1:]-r_obsdates[0:-1]

plt.figure(figsize=(6,4))
plt.hist(g_dt, range=(0,12), histtype='step', color='k', density=True)
plt.xlabel("Days before next obs", fontsize=16)
plt.ylabel("Fraction of $g$-band observations", fontsize=16)
plt.tick_params(labelsize=14)
plt.tight_layout()
plt.show()
