""" Estimate the rates of transients like AT2018cow and ZTF18abvkwla """

import numpy as np
import glob
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
from print_diagnostics import get_lc

# First, construct a template that is exactly its LC
jd,dt,mag,emag,filt,pid = get_lc('ZTF18abvkwla')
Mag = mag-Planck15.distmod(z=0.2714).value

# Luminosity function: ignore for now

# For each field, determine detection efficiency
# within 140 distance bins between z = 0.005 and z = 0.7
# for each bin, 
### choose a random distance within the bin
### k-correct template LC, assuming spec is blackbody
### choose random explosion epoch between XX and XX 
### map resulting explosion onto actual observation dates
### assign random noise pulled from actual distribution of noises...
### run resulting synthetic data through same LC selection criteria
### efficiency is the fraction of the 10,000 iterations that were recovered
