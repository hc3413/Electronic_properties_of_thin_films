
#import all the libraries needed
from import_dep import *

def set_plot_style(scale=1.0, use_tex=True):
    """
    Set publication-quality plot styles with scaling.
    
    Parameters:
    -----------
    scale : float
        Scaling factor for figure size and elements (default: 1.0)
    use_tex : bool
        Whether to use LaTeX for rendering text (default: True)
    """
    # Use a colorblind-friendly colormap with at least 10 distinct colors
    cmap_colors = sns.color_palette("colorblind", 10)
    
    # Science style settings with proper scaling
    plt.rcParams.update({
        # Figure settings
        'figure.figsize': [3.5 * scale, 2.625 * scale],
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05 * scale,
        
        # Font and text settings
        'font.family': ['serif'],
        'font.size': 9 * scale,  # Base font size
        'axes.labelsize': 10 * scale,
        'axes.titlesize': 10 * scale,
        'xtick.labelsize': 8 * scale,
        'ytick.labelsize': 8 * scale,
        'legend.fontsize': 8 * scale,
        'mathtext.fontset': 'dejavuserif',
        'text.usetex': use_tex,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{amssymb}',
        
        # Axes settings
        'axes.linewidth': 0.5 * scale,
        'axes.prop_cycle': cycler('color', ['#0C5DA5', '#00B945', '#FF9500', 
                                           '#FF2C00', '#845B97', '#474747', '#9e9e9e']),
        
        # Grid settings
        'grid.linewidth': 0.5 * scale,
        
        # Legend settings
        'legend.frameon': False,
        
        # Line settings
        'lines.linewidth': 1.0 * scale,
        'lines.markersize': 5.0 * scale,
        
        # Tick settings
        'xtick.direction': 'in',
        'xtick.major.size': 3.0 * scale,
        'xtick.major.width': 0.5 * scale,
        'xtick.minor.size': 1.5 * scale,
        'xtick.minor.visible': True,
        'xtick.minor.width': 0.5 * scale,
        'xtick.top': True,
        
        'ytick.direction': 'in',
        'ytick.major.size': 3.0 * scale,
        'ytick.major.width': 0.5 * scale,
        'ytick.minor.size': 1.5 * scale,
        'ytick.minor.visible': True,
        'ytick.minor.width': 0.5 * scale,
        'ytick.right': True,
        
    })
    
    return
