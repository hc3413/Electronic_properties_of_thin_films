#import all the libraries needed
from import_dep import *

def set_plot_style(fig_size, use_tex=True):
    """
    Set publication-quality plot styles.
    
    Parameters:
    -----------
    use_tex : bool
        Whether to use LaTeX for rendering text (default: True)
    """
    # Use a colorblind-friendly colormap with at least 10 distinct colors
    cmap_colors =   sns.color_palette("colorblind", 12) #sns.color_palette("bright", 10)
    color_cycler = cycler('color', cmap_colors)
    color_cycler_2 = cycler('color', ['#0C5DA5', '#00B945', '#FF9500', 
                                           '#FF2C00', '#845B97', '#474747', '#9e9e9e'])
    
    # Science style settings
    plt.rcParams.update({
        # Figure settings
        'figure.figsize':fig_size,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
        
        # Font and text settings
        'font.family': ['serif'],
        'font.size': 9,  # Base font size
        'axes.labelsize': 10,
        'axes.titlesize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'mathtext.fontset': 'dejavuserif',
        'text.usetex': use_tex,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{amssymb}',
        
        # Axes settings
        'axes.linewidth': 0.5,
        'axes.prop_cycle': color_cycler ,
        
        # Grid settings
        'grid.linewidth': 0.5,
        'axes.grid': True,
        
        # Legend settings
        'legend.frameon': False,
        
        # Line settings
        'lines.linewidth': 1.0,
        'lines.markersize': 4.0,
        # Errorbar settings
        'errorbar.capsize': 0,
        
        # Tick settings
        'xtick.direction': 'in',
        'xtick.major.size': 3.0,
        'xtick.major.width': 0.5,
        'xtick.minor.size': 1.5,
        'xtick.minor.visible': True,
        'xtick.minor.width': 0.5,
        'xtick.top': True,
        
        'ytick.direction': 'in',
        'ytick.major.size': 3.0,
        'ytick.major.width': 0.5,
        'ytick.minor.size': 1.5,
        'ytick.minor.visible': True,
        'ytick.minor.width': 0.5,
        'ytick.right': True,
        
        # Prevent autolayout to ensure that the figure size is obeyed
        'figure.autolayout': False,
    })
    
    return


def add_colorbar(fig, ax, sm, min_val, max_val, fig_size, field = True):
    """
    Add and adjust a colorbar to the given axis.
    
    Parameters:
    - fig: The figure object.
    - ax: The axis object to which the colorbar will be added.
    - sm: The ScalarMappable object for the colorbar.
    - min_field: The minimum value for the colorbar ticks.
    - max_field: The maximum value for the colorbar ticks.
    - fig_size: The size of the figure.
    - field: Whether the colorbar represents a magnetic field (True) or temperature (False).
    """
    cax = fig.add_subplot(ax)  # Use the provided axis for the colorbar
    cbar = plt.colorbar(sm, cax=cax)
    cbar.set_ticks([min_val, max_val])
    
    if field:
        cbar.set_ticklabels([f'{min_val:.1f} T', f'{max_val:.1f} T']) # Field
    else:
        cbar.set_ticklabels([f'{min_val:.1f} K', f'{max_val:.1f} K']) # Temperature
   
    cbar.minorticks_off()  # Remove minor ticks
    cbar.outline.set_linewidth(0.5)

    # Adjust colorbar position and size based on figure size
    if fig_size == [3.5, 2.625]:
        height_scale = 0.8
        pos = cax.get_position()
        new_height = pos.height * height_scale
        new_y0 = pos.y0 + (pos.height - new_height) / 2  # Center the colorbar vertically
        cax.set_position([
            pos.x0 + 0.03,  # Adjust x0 to move the colorbar to the right
            new_y0,  # Center the colorbar vertically
            pos.width * 0.3,  # Adjust width to shrink the colorbar
            new_height  # Adjust height to shrink the colorbar
        ])
    else:
        height_scale = 0.8
        pos = cax.get_position()
        new_height = pos.height * height_scale
        new_y0 = pos.y0 + (pos.height - new_height) / 2  # Center the colorbar vertically
        cax.set_position([
            pos.x0 + 0.03,  # Adjust x0 to move the colorbar to the right
            new_y0,  # Center the colorbar vertically
            pos.width * 0.3,  # Adjust width to shrink the colorbar
            new_height  # Adjust height to shrink the colorbar
        ])
