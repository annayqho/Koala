import matplotlib
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MultipleLocator,NullFormatter
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)

fig,ax = plt.subplots(1,1,figsize=(6,4))

# X-band
ax.errorbar(
        [81, 310, 352, 396], [0.364, 0.061, 0.045, 0.031], 
        yerr=[0.006, 0.003, 0.003, 0.003],
        c='k', marker='s', label="VLA: 10 GHz", ms=8)
ax.plot(
        [81, 310, 352, 396], [0.364, 0.061, 0.045, 0.031], 
        c='k', marker='s')

# C-band
ax.errorbar([343,397], [0.089,0.033], yerr=[0.003,0.003], 
        c='#f98e09', marker='D', ms=8,
        label="VLA: 6 GHz")
ax.plot([343,397], [0.089,0.033], c='#f98e09')

# S-band
ax.errorbar([346], [0.067], yerr=[0.005], mfc='white', 
        mec='k', marker='o', ms=9, c='k',
        label="VLA: 3 GHz")
ax.scatter([188], [0.134], facecolor='white', 
        edgecolor='k', marker='o', s=80, label="_none", zorder=5)
ax.plot([188,346], [0.134,0.067], c='k', ls='--')
ax.arrow(
        188, 0.134, 0, -0.03, color='k', 
        length_includes_head=True, head_width=10, head_length=0.01)

# L-band
ax.errorbar([351], [0.135], yerr=[0.007],
        c='#57106e', marker='^', ms=10,
        label="VLA: 1.5 GHz")

# GMRT point
ax.scatter(364, 0.105, c='k', marker='o', s=30, label="GMRT: 650 MHz")
ax.arrow(
        364, 0.105, 0, -0.02, color='k', 
        length_includes_head=True, head_width=15, head_length=0.01)


# zoom-in showing the SED
axins = inset_axes(
        ax, 2, 1, loc=3,
        bbox_to_anchor=(0.15,0.15),
        bbox_transform=ax.transAxes)
axins.errorbar(10, 0.045, 0.003, marker='s', c='k')
axins.errorbar(6, 0.089, 0.003, marker='D', c='#f98e09')
axins.errorbar(3, 0.067, 0.005, marker='o', mfc='white', mec='black')
axins.errorbar(1.5, 0.135, 0.007, marker='^', c='#57106e')
axins.scatter(0.65, 0.105, marker='o', c='k')
axins.arrow(0.65, 0.105, 0, -0.02, color='k', head_width=0.5, head_length=0.01)
axins.tick_params(axis='both', labelsize=12)
axins.set_ylabel("Flux (mJy)", fontsize=12)
axins.set_xlabel("$\\nu$ (GHz)", fontsize=12)
axins.text(0.9,0.9,r'$\Delta t\approx350$\,d', transform=axins.transAxes,
        fontsize=11, horizontalalignment='right', verticalalignment='top')

#axins.yaxis.set_label_position("right")
#axins.yaxis.tick_right()

ax.set_ylabel("Flux (mJy)", fontsize=12)
ax.set_xlabel("$\Delta t$ (d)", fontsize=12)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_yticks([0.1,0.3,0.05])
ax.set_xticks([100,200,300,400])
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.tick_params(axis='both', labelsize=12)
ax.legend(fontsize=10, loc='upper center', ncol=2)
ax.xaxis.set_minor_formatter(NullFormatter())

fig.tight_layout()
plt.savefig("radio_lc.png", dpi=300)
#plt.show()
