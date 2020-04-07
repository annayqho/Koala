""" Perform rate estimate with limiting magnitude cuts """

from astropy.io import ascii

# The fields that survived the cut last time
dat = ascii.read("../data/fieldnights_1DC.txt",delimiter=',')

# The limiting magnitudes of individual fields
maglims = np.loadtxt("seeing.txt", delimiter=',')

# Now, check which of these field-nights have a limiting magnitude of >19.75
choose = maglims[:,1] > 19.75

# What field-nights does this correspond to?
ind = maglims[:,0][choose].astype(int)
fn_shallow = dat[ind]

choose = maglims[:,1] > 20.5

# What field-nights does this correspond to?
ind = maglims[:,0][choose].astype(int)
fn_deep = dat[ind]
