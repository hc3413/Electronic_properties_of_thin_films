#import all the libraries needed
from import_dep import *


def set_plot_style(export_data = False, powerpoint_data = False, use_tex=True):
    """
    Set publication-quality plot styles.
    
    Parameters:
    -----------
    use_tex : bool
        Whether to use LaTeX for rendering text (default: True)
    """
    
    # Set the figure size based on whether we are visualising or exporting the data
    if export_data == True:
        fig_size = [3.5, 2.625] # Publication ready sizes
    else:
        fig_size = [9, 6] # Better for visualisation
    
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
        'legend.fontsize': 7,
        'mathtext.fontset': 'dejavuserif',
        'text.usetex': use_tex,
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{amssymb}',
        
        # Axes settings
        'axes.linewidth': 0.5,
        'axes.prop_cycle': color_cycler ,
        
        # Grid settings
        'grid.linewidth': 0.5,
        'axes.grid': True,
        'axes.axisbelow': True,
        
        # Legend settings
        'legend.frameon': True,
        'legend.framealpha': 0.4,
        
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
    
    return fig_size


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


from matplotlib.colors import Normalize
import numpy as np
import matplotlib.pyplot as plt

def generate_colormaps_and_normalizers(dat):
    """
    Generate colormaps and normalizers for temperature and field values.
    Taking the input data to work out the maximum and minimum values for the colormaps.
    Parameters:
    -----------
    dat : list
        List of data objects containing temperature and field information.

    Returns:
    --------
    cmap_temp : Colormap
        Colormap for temperature values.
    cmap_field : Colormap
        Colormap for field values.
    norm_temp : Normalize
        Normalizer for temperature values.
    norm_field : Normalize
        Normalizer for field values.
    cmap_dat : ndarray
        Colormap for distinguishing between datasets.
    """
    # Extract raw temperature values to prevent rounding
    all_temps = []
    for d in dat:
        temps = np.copy(d.tf_av).reshape((d.ctf[4] * d.ctf[5], 2))
        all_temps = np.concatenate([all_temps, temps[:, 0]])  # Extract the temperature values (first column)

    # Concatenate all field arrays
    all_fields = np.concatenate([d.ctf[2] for d in dat])

    # Find the min and max values
    min_temp = np.min(all_temps)
    max_temp = np.max(all_temps)
    min_field = np.min(all_fields)
    max_field = np.max(all_fields)

    # Normalize the temperature and field values
    norm_temp = Normalize(vmin=min_temp, vmax=max_temp)
    norm_field = Normalize(vmin=min_field, vmax=max_field)

    # Generate colormaps
    cmap_temp = plt.get_cmap('coolwarm')
    cmap_field = plt.get_cmap('coolwarm')
    
    # Generate a list of markers for the data
    mark_p = [ 'x', 'o', '*', 'd', '^', 'v','+',  '<', '>', 'p', 'P', 'h', 'H', 'X', 'D', '|', '_', '1', '2', '3', '4', '8', 's', 'p', 'P', 'o', 'h', 'H', 'X', 'd', 'D', '|', '_', '1', '2', '3', '4', '8', 's']

   

    return cmap_temp, cmap_field, norm_temp, norm_field, mark_p, min_temp, max_temp, min_field, max_field




# Define a function for the PPTX module to add a slide with a title and image
# --- Define the standard font for the slide (consistent usage) ---
SLIDE_FONT = "Avenir" # Note: Ensure this font is available on the system running the code AND the system viewing the PPTX

def add_slide(fig, title, notes, prs, path_out):
    """
    Adds a slide to the presentation with a custom title bar (gradient),
    figure on the left, and notes on the right, using the defined SLIDE_FONT.

    Args:
        fig: The matplotlib figure object to add.
        title (str): The title for the slide.
        notes (list): A list of strings for the notes section.
        prs: The Presentation object.
        path_out (str or Path): The directory to save the temporary image file.
    """
    # Check if figure exists and is valid
    if fig is None or not hasattr(fig, 'axes') or len(fig.axes) == 0:
        print(f"Skipping '{title}': Figure does not exist or is empty.")
        if fig:
            try:
                import matplotlib.pyplot as plt
                plt.close(fig)
            except ImportError: pass
            except Exception as e_close: print(f"Error closing invalid figure '{title}': {e_close}")
        return

    # --- Use BLANK layout ---
    slide_layout = prs.slide_layouts[6]
    slide = prs.slides.add_slide(slide_layout)

    # --- Set background color ---
    background = slide.background
    fill = background.fill
    fill.solid()
    fill.fore_color.rgb = RGBColor(240, 240, 240)

    # --- Create Title Rectangle ---
    slide_width = prs.slide_width
    slide_height = prs.slide_height
    title_box_height = Inches(0.819)
    title_box = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE,
        left=0, top=0, width=slide_width, height=title_box_height
    )
    # Remove border
    title_box.line.fill.background()
    title_box.line.width = Pt(0)

    # --- Apply Linear Gradient Fill to Title Rectangle ---
    fill = title_box.fill
    fill.gradient()
    fill.gradient_angle = 0.0 # Linear Left to Right

    # --- Define gradient stops ---
    if len(fill.gradient_stops) >= 2: # Check if stops exist
        stop1 = fill.gradient_stops[0]
        stop1.color.rgb = RGBColor(40, 80, 160) # Left color
        stop1.position = 0.1

        stop2 = fill.gradient_stops[1]
        stop2.color.rgb = RGBColor(100, 160, 220) # Right color
        stop2.position = 1.0
    else:
        print(f"Warning: Could not set gradient stops for title '{title}'. Default gradient may apply.")
        # Fallback: fill.solid(); fill.fore_color.rgb = RGBColor(68, 114, 196)

    # --- Add Text Box FOR the title ON TOP ---
    title_left = Inches(0.2)
    title_top = Inches(0.1)
    title_width = slide_width - Inches(0.4)
    title_textbox_height = title_box_height - Inches(0.1)
    title_textbox = slide.shapes.add_textbox(
        title_left, title_top, title_width, title_textbox_height
    )
    tf = title_textbox.text_frame
    tf.margin_left = 0
    tf.margin_right = 0
    tf.margin_top = 0
    tf.margin_bottom = 0
    tf.word_wrap = False
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE

    p = tf.paragraphs[0]
    p.text = title
    p.font.name = SLIDE_FONT
    p.font.size = Pt(20)
    p.font.bold = False
    p.font.color.rgb = RGBColor(255, 255, 255)
    p.alignment = PP_ALIGN.LEFT

    # --- Save the figure as an image ---
    output_dir = Path(path_out)
    output_dir.mkdir(parents=True, exist_ok=True)
    safe_title = "".join(c if c.isalnum() else "_" for c in title)
    img_path = output_dir / f'{safe_title}.png'

    try:
        fig.savefig(img_path, dpi=300, bbox_inches='tight', transparent=True)
        try:
            import matplotlib.pyplot as plt
            plt.close(fig)
        except ImportError: print("Warning: Matplotlib not imported, cannot close figure.")
        except Exception as e_close: print(f"Error closing figure '{title}': {e_close}")
    except Exception as e_save:
        print(f"Error saving figure '{title}': {e_save}")
        if fig:
            try:
                import matplotlib.pyplot as plt
                plt.close(fig)
            except Exception: pass
        return

    # --- Define Layout Areas ---
    margin = Inches(0.2)
    content_top = title_box_height + margin
    content_height = slide_height - content_top - margin
    # Calculate total width available for columns after accounting for 3 margins (left, center, right)
    total_content_width = slide_width - (3 * margin)
    # Left column takes 60% of the available width
    left_column_width = 0.6 * total_content_width
    # Right column takes 40% of the available width
    right_column_width = 0.4 * total_content_width
    # Right column starts after the left margin, the left column width, and the center margin
    right_column_left = margin + left_column_width + margin

    # --- Add the image to the left column ---
    img_left = margin
    img_top = content_top
    max_img_width = left_column_width
    max_img_height = content_height

    try:
        # Add picture, initially using max width
        # Ensure initial position and width are integers (though Inches usually handles this)
        pic = slide.shapes.add_picture(str(img_path), int(img_left), int(img_top), width=int(max_img_width))

        # Scale height if it exceeds max height, adjusting width proportionally
        if pic.height > max_img_height:
            scale_factor = max_img_height / pic.height
            pic.height = int(max_img_height) # Cast to int
            pic.width = int(pic.width * scale_factor) # Cast to int
            # Recenter horizontally after scaling
            pic.left = int(img_left + (max_img_width - pic.width) / 2) # Cast to int

        # Scale width if it exceeds max width (less likely after initial width setting, but good practice)
        elif pic.width > max_img_width: # Use elif to avoid double scaling if both are too large initially
             scale_factor = max_img_width / pic.width
             pic.width = int(max_img_width) # Cast to int
             pic.height = int(pic.height * scale_factor) # Cast to int
             # Recenter vertically after scaling
             pic.top = int(img_top + (max_img_height - pic.height) / 2) # Cast to int

        # # Center vertically if image is smaller than available height
        # elif pic.height < max_img_height:
        #      pic.top = int(img_top + (max_img_height - pic.height) / 2) # Cast to int
        
        # # Center horizontally if image is smaller than available width (added for completeness)
        # elif pic.width < max_img_width:
        #      pic.left = int(img_left + (max_img_width - pic.width) / 2) # Cast to int


        print(f"Added slide: {title}")

        # --- Add Notes Box to the right column if notes are provided ---
        if notes and isinstance(notes, list) and any(n and n.strip() for n in notes):
            notes_left = right_column_left
            notes_top = content_top
            notes_width = right_column_width
            notes_height = content_height # Use full available height

            # Ensure notes box dimensions are integers
            notes_box = slide.shapes.add_textbox(int(notes_left), int(notes_top), int(notes_width), int(notes_height))
            notes_tf = notes_box.text_frame
            notes_tf.word_wrap = True
            notes_tf.vertical_anchor = MSO_ANCHOR.TOP # Anchor text to top
            notes_tf.clear()

            for line in notes:
                line = line.strip()
                if line:
                    p_notes = notes_tf.add_paragraph()
                    # Add bullet point using unicode character
                    p_notes.text = f"• {line}"
                    p_notes.font.name = SLIDE_FONT
                    p_notes.font.size = Pt(14)
                    p_notes.font.color.rgb = RGBColor(0, 0, 0)
                    p_notes.line_spacing = 1.1
                    p_notes.space_before = Pt(3) # Add space before bullet point
                    p_notes.space_after = Pt(3) # Add space after bullet point

            # Remove extra space after the last paragraph
            if len(notes_tf.paragraphs) > 0:
                notes_tf.paragraphs[-1].space_after = Pt(0)

    except FileNotFoundError:
        print(f"Error adding picture: Image file not found at '{img_path}'.")
    except Exception as e_add:
        print(f"Error adding picture or notes for slide '{title}': {e_add}")

    # Clean up the saved image file
    try:
        os.remove(img_path)
    except OSError as e_remove:
        print(f"Warning: Could not remove temporary image file '{img_path}': {e_remove}")


def sample_tag_gen(dat):
    '''function to generate a list of sample tags for each sample in the dat list 
    in the format [dat[0].sample_code = dat[0].material, dat[1].sample_code = dat[1].material, ...]'''
    tag_list = []
    for d in dat:
        tag_list.append(f"{d.sample_code} = {d.material[1:-1]}")
    return tag_list
