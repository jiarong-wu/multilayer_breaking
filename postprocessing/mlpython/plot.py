import matplotlib as mpl
from matplotlib import pyplot as plt

def set_letters(x=-0.2, y=1.05, fontsize=11, letters=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p'], color='k'):
    fig = plt.gcf()
    axes = fig.axes
    j = 0
    for ax in axes:
        if hasattr(ax, 'collections'):
            if len(ax.collections) > 0:
                collection = ax.collections[0]
            else:
                collection = ax.collections
            if isinstance(collection, mpl.collections.LineCollection):
                print('Colorbar-like object skipped')
            else:
                try:
                    ax.text(x,y,f'({letters[j]})', transform = ax.transAxes, fontweight='bold', fontsize=fontsize, color=color)
                except:
                    print('Cannot set letter', letters[j])
                j += 1
                
def add_colorbar (fig, loc=[0.15,0.62,0.01,0.15], cmap='Oranges', vmax=0.2, vmin=0.05, ticks=(0.05,0.1,0.15,0.2)):
    ax2 = fig.add_axes(loc)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cbar = mpl.colorbar.ColorbarBase(ax2, cmap=plt.get_cmap(cmap), norm=norm, orientation='vertical', ticks=ticks)
    cbar.ax.tick_params(labelsize=6)
    # cbar.ax.text(0.5, 1.15, r'$\sigma$', ha='center', va='center', transform=cbar.ax.transAxes)
    return cbar