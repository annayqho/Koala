import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)

# four frequencies
fig,axarr = plt.subplots(2,2,figsize=(5,5))

# X-band
ax = axarr[0]
ax.scatter([81, 310, 352, 396], [0.364, 0.061, 0.045, 0.031], c='k', marker='o')
ax.text(0,0,'10 GHz',fontsize=10,transform=ax.transAxes)

# C-band
ax = axarr[0]
ax.scatter([81, 310, 352, 396], [0.364, 0.061, 0.045, 0.031], c='k', marker='o')

fig.tight_layout()
plt.show()
