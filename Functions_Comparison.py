# Import all the libraries needed
from import_dep import *
from Class_Import import *
from Functions_style import set_plot_style

from typing import Optional, List, Tuple, Union, Dict
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches


# =============================================================================
# Material definitions (shared with RSM_functions.py)
# =============================================================================
PEROVSKITE_MATERIALS = {
    'STO':      {'name': 'SrTiO3',              'a_bulk': 3.905},
    'BSO':      {'name': 'BaSnO3',              'a_bulk': 4.116},
    'SSO':      {'name': 'SrSnO3',              'a_bulk': 4.036},
    'LSO':      {'name': 'LaScO3',              'a_bulk': 4.050},
    'BTO':      {'name': 'BaTiO3',              'a_bulk': 4.010},
    'BSSO_25':  {'name': 'Ba0.25Sr0.75SnO3',    'a_bulk': 4.056},
    'BSSO_40':  {'name': 'Ba0.40Sr0.60SnO3',    'a_bulk': 4.068},
    'BSSO_55':  {'name': 'Ba0.55Sr0.45SnO3',    'a_bulk': 4.080},
    'BSSO_70':  {'name': 'Ba0.70Sr0.30SnO3',    'a_bulk': 4.092},
}

MATERIAL_MARKERS = {
    'STO': '*', 'BSO': 'o', 'SSO': 's', 'LSO': '^', 'BTO': 'd',
    'BSSO': 'X',
}


# =============================================================================
# Helper functions
# =============================================================================
def _get_material_marker(mat: str, default: str = 'x') -> str:
    """Return marker for a material, matching BSSO prefix for all BSSO_xx variants."""
    if mat in MATERIAL_MARKERS:
        return MATERIAL_MARKERS[mat]
    for prefix in MATERIAL_MARKERS:
        if mat.startswith(prefix + '_'):
            return MATERIAL_MARKERS[prefix]
    return default


def _material_matches(mat_name: str, select_list: List[str]) -> bool:
    """Check if *mat_name* matches any entry in *select_list* (prefix-aware)."""
    mat_upper = mat_name.upper()
    for s in select_list:
        s_upper = s.upper()
        if mat_upper == s_upper or mat_upper.startswith(s_upper + '_'):
            return True
    return False


def _resolve_sample_codes(sample_list: list) -> List[str]:
    """
    Resolve a mixed list of PPMSData objects and/or sample code strings
    to a list of uppercase sample code strings.
    """
    codes = []
    for item in sample_list:
        if hasattr(item, 'sample_code'):
            codes.append(item.sample_code.upper())
        elif isinstance(item, str):
            codes.append(item.upper())
        else:
            codes.append(str(item).upper())
    return codes


def _get_ppms_object(sample_list: list, code: str):
    """Return the PPMSData object matching a sample code, or None."""
    code_upper = code.upper()
    for item in sample_list:
        if hasattr(item, 'sample_code') and item.sample_code.upper() == code_upper:
            return item
    return None


# =============================================================================
# Database loading
# =============================================================================
def load_growth_database(
    xlsx_path: str,
    sample_codes: Optional[List[str]] = None,
    layer_material: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Load growth condition and XRD data from the shared Excel peak database.

    Each sample has multiple rows (one per layer: substrate, BSO, SSO, LSO, etc.).
    This function retains the per-layer structure so that growth conditions and
    XRD parameters specific to each layer are preserved.

    Parameters:
        xlsx_path      : Path to the .xlsx file (e.g., RSM_peaks.xlsx).
        sample_codes   : Optional list of sample codes to filter by (case-insensitive).
        layer_material : Optional list of material names to filter layers by
                         (e.g., ['BSO'] to keep only the BSO layer for each sample).
                         Uses prefix matching: 'BSSO' matches 'BSSO_25', 'BSSO_40', etc.
                         If None, all non-substrate layers are returned.

    Returns:
        DataFrame with one row per (sample_code, material layer), containing
        both growth conditions and XRD parameters for that specific layer.
    """
    xlsx_path = Path(xlsx_path)
    if not xlsx_path.exists():
        print(f"File not found: {xlsx_path}")
        return pd.DataFrame()

    df = pd.read_excel(xlsx_path, engine='openpyxl')

    # Filter by sample codes if provided
    if sample_codes is not None:
        codes_upper = [c.upper() for c in sample_codes]
        if 'sample_code' in df.columns:
            df = df[df['sample_code'].str.upper().isin(codes_upper)]
        else:
            print("Warning: 'sample_code' column not found in database.")
            return pd.DataFrame()

    if df.empty:
        return df

    # Remove substrate rows
    if 'is_substrate' in df.columns:
        df = df[df['is_substrate'] != True].copy()

    # Filter by layer material if specified
    if layer_material is not None and 'material' in df.columns:
        df = df[df['material'].apply(
            lambda m: _material_matches(m, layer_material) if pd.notna(m) else False
        )].copy()

    # Keep relevant columns
    growth_cols = [
        'sample_code', 'material', 'Ba_doping', 'La_doping',
        'O2_growth_pressure', 'T_growth', 'n_pulses', 'thickness_nm', 'Hz', 'E_mJ',
    ]
    xrd_cols = ['a', 'c', 'qx', 'qz', 'intensity', 'intensity_measured',
                'intensity_norm', 'sigma_x', 'sigma_z', 'eta', 'background']

    available_cols = [c for c in growth_cols + xrd_cols if c in df.columns]
    df = df[available_cols].copy()

    # Compute derived columns
    if 'a' in df.columns and 'c' in df.columns:
        df['tetragonality'] = df['c'] / df['a']
    if 'intensity_norm' in df.columns and 'n_pulses' in df.columns:
        df['intensity_norm_perpulse'] = df['intensity_norm'] / df['n_pulses']

    return df


# =============================================================================
# Prompt for missing data
# =============================================================================
def _prompt_for_missing(df: pd.DataFrame, x_col: str, x_axis_var: str,
                        layer_material: List[str],
                        sample_list: Optional[list] = None) -> pd.DataFrame:
    """
    Check for samples with missing x-axis values and prompt the user to input them.

    Only prompts once per unique sample_code (the entered value is applied to all
    temperature rows for that sample).  Entered values are also stored back into
    the PPMSData object's ``sample_parameters`` dict so they persist across calls.

    Parameters:
        df             : Merged DataFrame with electrical + growth data.
        x_col          : The actual column name used for the x-axis.
        x_axis_var     : The user-facing name of the x-axis variable.
        layer_material : The material layer(s) selected for growth conditions.
        sample_list    : Optional list of PPMSData objects — used to persist
                         prompted values into ``sample_parameters``.

    Returns:
        DataFrame with missing values filled in from user input.
    """
    missing_mask = df[x_col].isna()
    if not missing_mask.any():
        return df

    missing_codes = df.loc[missing_mask, 'sample_code'].unique()

    layer_str = ', '.join(layer_material) if layer_material else 'all layers'
    print(f"\n--- Missing '{x_axis_var}' values for {len(missing_codes)} sample(s) "
          f"(layer: {layer_str}) ---")

    df = df.copy()
    for code in missing_codes:
        prompt_str = (f"  Enter '{x_axis_var}' for {code} ({layer_str}) "
                      f"[press Enter to skip]: ")
        user_input = input(prompt_str)
        if user_input.strip() == '':
            print(f"    Skipped {code} — will be excluded from plot.")
            continue
        try:
            value = float(user_input.strip())
            df.loc[df['sample_code'] == code, x_col] = value
            print(f"    Set {code} '{x_axis_var}' = {value}")

            # Persist the entered value into the PPMSData object
            if sample_list is not None:
                ppms = _get_ppms_object(sample_list, code)
                if ppms is not None:
                    for layer in (layer_material if layer_material else ['Unknown']):
                        if layer not in ppms.sample_parameters:
                            ppms.sample_parameters[layer] = {}
                        ppms.sample_parameters[layer][x_axis_var] = value
        except ValueError:
            print(f"    Invalid input '{user_input}' — {code} will be excluded from plot.")

    return df


# =============================================================================
# Extract electrical data at specific temperatures
# =============================================================================
def _extract_electrical_data(
    sample_list: list,
    temperatures: List[float],
    temp_tolerance: float = 15.0,
) -> pd.DataFrame:
    """
    Extract electrical transport data from PPMSData objects at specified temperatures.

    For each sample and each requested temperature, finds the closest measured
    temperature within `temp_tolerance` and extracts:
      - resistivity (zero-field or averaged), sheet resistance
      - mobility, carrier density, Hall coefficient

    Parameters:
        sample_list    : List of PPMSData objects.
        temperatures   : List of temperatures (K) to extract data at.
        temp_tolerance : Maximum allowed difference (K) between requested and
                         measured temperature (default: 5 K).

    Returns:
        DataFrame with columns:
            sample_code, temperature, rho_xx, rho_xx_error, sheet_resistance,
            mobility, mobility_error, carrier_density, carrier_density_error,
            hall_coefficient, material
    """
    rows = []

    for ppms in sample_list:
        if not hasattr(ppms, 'sample_code') or ppms.sample_code is None:
            continue

        code = ppms.sample_code.upper()
        mat = ppms.material if ppms.material else 'Unknown'
        thickness = ppms.film_thickness

        for T_req in temperatures:
            row = {
                'sample_code': code,
                'temperature': T_req,
                'material': mat,
            }

            # Track the actual measured temperature (closest match)
            actual_temp = T_req
            _hall_found = False
            _res_found = False

            # --- Extract from hall_coefficient array ---
            # VDP columns: 'Temp (K)', Hallco_A, R^2_A, Hallco_B, R^2_B, Hallco_avg, R^2_avg,
            #              n (cm^-2 or cm^-3), n_err, mobility, mobility_err, n_fitted, mu_fitted
            # HallBar columns: same structure but different header names
            if ppms.hall_coefficient is not None and ppms.hall_coefficient.shape[0] > 0:
                hc = ppms.hall_coefficient
                temps_measured = hc[:, 0]
                idx = np.argmin(np.abs(temps_measured - T_req))
                diff = np.abs(temps_measured[idx] - T_req)
                if diff <= temp_tolerance:
                    actual_temp = float(temps_measured[idx])
                    row['hall_coefficient'] = hc[idx, 5]   # average Hall coefficient
                    row['carrier_density'] = hc[idx, 7]    # n
                    row['carrier_density_error'] = hc[idx, 8]  # n_error
                    row['mobility'] = hc[idx, 9]           # mobility (cm^2/Vs)
                    row['mobility_error'] = hc[idx, 10]    # mobility_error
                    _hall_found = True
                    print(f"  {code}: Hall data at {actual_temp:.1f} K "
                          f"(requested {T_req} K, Δ={diff:.1f} K)")
                else:
                    print(f"  {code}: Hall — closest T is {temps_measured[idx]:.1f} K "
                          f"(requested {T_req} K, Δ={diff:.1f} K > tolerance {temp_tolerance} K) — SKIPPED")
            else:
                print(f"  {code}: No hall_coefficient data available")

            # --- Extract from res_data array ---
            # Columns: Temp(K), Field(T), rho_A, rho_B, rho_avg, rho_err, rho_fit
            if ppms.res_data is not None and ppms.res_data.shape[0] > 0:
                rd = ppms.res_data  # shape: (n_temp, n_field, 7)
                if rd.ndim == 3:
                    # Find the temperature index closest to T_req
                    temps_res = rd[:, 0, 0]  # temperatures at first field point
                    idx_t = np.argmin(np.abs(temps_res - T_req))
                    diff_res = np.abs(temps_res[idx_t] - T_req)
                    if diff_res <= temp_tolerance:
                        # Find zero-field (or closest to zero) data point
                        fields = rd[idx_t, :, 1]
                        idx_f = np.argmin(np.abs(fields))
                        rho = rd[idx_t, idx_f, 4]    # rho_xx_average
                        rho_err = rd[idx_t, idx_f, 5]  # rho_error
                        row['rho_xx'] = rho
                        row['rho_xx_error'] = rho_err
                        # Sheet resistance = rho / thickness
                        if thickness is not None and thickness != 1:
                            row['sheet_resistance'] = rho / thickness
                        else:
                            # If thickness == 1, res_data already IS sheet resistance
                            row['sheet_resistance'] = rho
                        _res_found = True
                        print(f"  {code}: Resistivity at {temps_res[idx_t]:.1f} K, "
                              f"B={fields[idx_f]:.2f} T "
                              f"(requested {T_req} K, Δ={diff_res:.1f} K)")
                    else:
                        print(f"  {code}: Resistivity — closest T is {temps_res[idx_t]:.1f} K "
                              f"(requested {T_req} K, Δ={diff_res:.1f} K > tolerance {temp_tolerance} K) — SKIPPED")
            else:
                print(f"  {code}: No res_data available")

            if not _hall_found and not _res_found:
                print(f"  {code}: ** No data extracted at {T_req} K — row will be empty **")

            # Use the REQUESTED temperature for grouping (colour/line/legend).
            # Points matched within tolerance all represent the same requested T.
            row['temperature'] = round(T_req)

            rows.append(row)

    return pd.DataFrame(rows)


# =============================================================================
# Main comparison plot function
# =============================================================================
def growth_optimisation_plot(
    sample_list: list,
    database_path: str,
    temperatures: List[float] = [300],
    layer_material: List[str] = ['BSO'],
    x_axis_var: str = 'O2_growth_pressure',
    y_axis_var: str = 'mobility',
    y_axis_var_right: Optional[str] = None,
    temp_tolerance: float = 5.0,
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    y_lim_right: Optional[Tuple[float, float]] = None,
    export_data: bool = False,
    fig_format: str = 'tiff',
    show_bulk_refs: bool = False,
    material_select: Union[bool, List[str]] = False,
    materials: Optional[Dict] = None,
    show_key: Union[bool, str] = True,
    line_plot: bool = False,
    fit_line: bool = False,
    log_y: bool = False,
    log_y_right: bool = False,
    log_x: bool = False,
) -> plt.Figure:
    """
    Growth-parameter optimisation plot combining growth conditions from the
    XRD Excel database with electrical transport data from PPMSData objects.

    Plots a chosen independent variable (growth condition) on the x-axis against
    a dependent variable (electrical property, XRD parameter, or growth condition)
    on the y-axis, with data points coloured/marked by temperature.

    Independent variables (x_axis_var):
        Growth conditions: 'Ba_doping', 'La_doping', 'O2_growth_pressure',
            'T_growth', 'n_pulses', 'thickness_nm', 'Hz', 'E_mJ'
        Sample identifier: 'sample_name'
        Manual input: 'input_new' — prompts for an axis title then a
            numeric value for each sample (one-off custom x-axis).
        (or any numeric column in the database)

    Dependent variables (y_axis_var / y_axis_var_right):
        Electrical:
            'mobility'          - Hall mobility (cm²/Vs)
            'carrier_density'   - Charge carrier density (cm⁻² or cm⁻³)
            'rho_xx'            - Volume resistivity (Ω·m)
            'sheet_resistance'  - Sheet resistance (Ω/□)
            'hall_coefficient'  - Hall coefficient
        XRD / structural:
            'a'                 - In-plane lattice parameter (Å)
            'c'                 - Out-of-plane lattice parameter (Å)
            'a_and_c'           - Both a and c on same axis
            'tetragonality'     - c/a ratio
            'intensity_norm'    - Normalised XRD intensity
            'intensity_norm_perpulse' - Normalised intensity per pulse
            'sigma_x', 'sigma_z', 'eta', 'background'
        Growth conditions (for cross-comparison):
            'T_growth', 'n_pulses', 'thickness_nm', etc.

    Parameters:
        sample_list      : PPMSData objects (with hall/resistivity data computed).
        database_path    : Path to the .xlsx peak/growth-condition database.
        temperatures     : List of temperatures (K) at which to extract electrical data.
                           Each temperature produces a differently coloured set of points.
        layer_material   : List of material layer names whose growth conditions and
                           XRD data to use from the database (default: ['BSO']).
                           Prefix-aware: ['BSSO'] matches BSSO_25, BSSO_40, etc.
                           The x-axis value and any XRD y-axis value are taken from
                           the row matching this material for each sample.
        x_axis_var       : Column name for x-axis.
        y_axis_var       : Column name for left y-axis.
        y_axis_var_right : Column name for right y-axis (optional).
        temp_tolerance   : Max T difference (K) for matching (default: 5).
        x_lim, y_lim, y_lim_right : Axis limits.
        export_data      : Use publication-ready figure sizing.
        fig_format       : Export format ('tiff', 'png', 'pdf', ...).
        show_bulk_refs   : Show horizontal lines for bulk lattice parameters.
        material_select  : False to plot all, or list of material labels to include.
        materials        : Material dictionary (default: PEROVSKITE_MATERIALS).
        show_key         : Legend control (True, False, or location string).
        line_plot        : Connect points with lines (sorted by x).
        fit_line         : Overlay linear regression fit line.
        log_y            : Use log scale for left y-axis.
        log_y_right      : Use log scale for right y-axis.
        log_x            : Use log scale for x-axis.

    Returns:
        matplotlib Figure.
    """
    if materials is None:
        materials = PEROVSKITE_MATERIALS

    # Resolve sample codes from mixed list of PPMSData objects and strings
    sample_codes = _resolve_sample_codes(sample_list)
    ppms_objects = [item for item in sample_list if hasattr(item, 'sample_code')]

    # ------------------------------------------------------------------
    # Build growth DataFrame: prefer cached sample_parameters, fall back
    # to the Excel database, then populate the cache for future calls.
    # ------------------------------------------------------------------
    cached_rows = []
    codes_needing_db = []

    for ppms in ppms_objects:
        code = ppms.sample_code.upper()
        has_cached = False
        if hasattr(ppms, 'sample_parameters') and ppms.sample_parameters:
            for layer, params in ppms.sample_parameters.items():
                if layer_material is None or _material_matches(layer, layer_material):
                    row = {'sample_code': code, 'material': layer}
                    row.update(params)
                    cached_rows.append(row)
                    has_cached = True
        if not has_cached:
            codes_needing_db.append(code)

    # Load from Excel only for samples not yet cached
    df_growth = pd.DataFrame()
    if codes_needing_db:
        # Load ALL layers from the database so every layer is cached for
        # future calls with a different layer_material selection.
        df_db_all = load_growth_database(database_path, sample_codes=codes_needing_db,
                                         layer_material=None)

        # Populate sample_parameters on each PPMSData object from the database
        if not df_db_all.empty:
            for ppms in ppms_objects:
                code = ppms.sample_code.upper()
                rows_for_sample = df_db_all[df_db_all['sample_code'].str.upper() == code]
                if rows_for_sample.empty:
                    continue
                if not hasattr(ppms, 'sample_parameters') or ppms.sample_parameters is None:
                    ppms.sample_parameters = {}
                for _, db_row in rows_for_sample.iterrows():
                    layer = db_row.get('material', 'Unknown')
                    if layer not in ppms.sample_parameters:
                        ppms.sample_parameters[layer] = {}
                    for col_name, col_val in db_row.items():
                        if col_name not in ('sample_code', 'material') and pd.notna(col_val):
                            ppms.sample_parameters[layer][col_name] = col_val

        # Filter to the requested layer_material for plotting
        if layer_material is not None and not df_db_all.empty and 'material' in df_db_all.columns:
            df_db = df_db_all[df_db_all['material'].apply(
                lambda m: _material_matches(m, layer_material) if pd.notna(m) else False
            )].copy()
        else:
            df_db = df_db_all

        df_growth = df_db
    elif not cached_rows:
        # All codes are strings (not PPMSData objects) — load everything from DB
        df_growth = load_growth_database(database_path, sample_codes=sample_codes,
                                         layer_material=layer_material)

    # Combine cached and freshly-loaded data
    if cached_rows:
        df_cached = pd.DataFrame(cached_rows)
        df_growth = pd.concat([df_cached, df_growth], ignore_index=True)

    # Extract electrical data at each requested temperature
    df_elec = _extract_electrical_data(ppms_objects, temperatures, temp_tolerance)

    if df_elec.empty:
        print("No electrical data found for the specified samples.")
        return plt.figure()

    # Merge: growth conditions + electrical data on sample_code
    # Each sample should now have one growth-data row (the selected layer),
    # replicated for each temperature from the electrical data.
    # Samples missing from the database will have NaN growth columns and
    # be handled by the missing-data prompt below.
    df_elec['sample_code'] = df_elec['sample_code'].str.upper()
    if not df_growth.empty:
        df_growth['sample_code'] = df_growth['sample_code'].str.upper()
        df_growth = df_growth.rename(columns={'material': 'layer_material'})
        df = pd.merge(df_elec, df_growth, on='sample_code', how='left')
    else:
        df = df_elec.copy()

    if df.empty:
        print("No data after merging electrical and growth data.")
        return plt.figure()

    # --- Handle special x-axis modes ---
    _sample_name_mode = (x_axis_var == 'sample_name')
    _input_new_mode = (x_axis_var == 'input_new')
    _input_new_label = None  # will hold user-supplied axis title

    if _input_new_mode:
        # Prompt for axis title
        _input_new_label = input("Enter x-axis title (e.g. 'time (hr)'): ").strip()
        if not _input_new_label:
            _input_new_label = 'Custom x'
        x_col = '_input_new_x'
        # Collect x-values in a plain list — one entry per PPMSData object,
        # correlated by list position.  Values are NOT persisted anywhere.
        _input_x_values = []  # same length as ppms_objects; None = skip
        for ppms in ppms_objects:
            code = ppms.sample_code.upper()
            tag = code
            if ppms.filename:
                tag += f' ({ppms.filename})'
            elif ppms.plot_str and ppms.plot_str != code:
                tag += f' ({ppms.plot_str})'
            prompt_str = (f"  Enter '{_input_new_label}' value for {tag} "
                          f"[press Enter to skip]: ")
            user_input = input(prompt_str)
            if user_input.strip() == '':
                print(f"    Skipped {tag} — will be excluded from plot.")
                _input_x_values.append(None)
                continue
            try:
                value = float(user_input.strip())
                _input_x_values.append(value)
                print(f"    Set {tag} '{_input_new_label}' = {value}")
            except ValueError:
                print(f"    Invalid input '{user_input}' — {tag} will be excluded.")
                _input_x_values.append(None)
        # Build df directly: extract electrical data per object and tag with
        # the x-value from the list.  Bypasses the growth merge entirely so
        # multiple objects sharing a sample_code each get their own rows.
        dfs = []
        for i, ppms in enumerate(ppms_objects):
            if _input_x_values[i] is None:
                continue
            df_one = _extract_electrical_data([ppms], temperatures, temp_tolerance)
            df_one[x_col] = _input_x_values[i]
            dfs.append(df_one)
        if not dfs:
            print('No values entered — nothing to plot.')
            return plt.figure()
        df = pd.concat(dfs, ignore_index=True)
    elif _sample_name_mode:
        _sample_order = sample_codes
        _x_map = {code: i for i, code in enumerate(_sample_order)}
        df = df.copy()
        df['_x_pos'] = df['sample_code'].str.upper().map(_x_map)
        df = df.dropna(subset=['_x_pos'])
        x_col = '_x_pos'
    else:
        x_col = x_axis_var
        if x_col not in df.columns:
            # Column doesn't exist at all — create it and prompt for all samples
            df[x_col] = np.nan
        # Prompt for missing x-axis values
        df = _prompt_for_missing(df, x_col, x_axis_var, layer_material,
                                 sample_list=sample_list)

    # Filter by material_select
    if material_select is not False and material_select:
        df = df[
            df['material'].apply(
                lambda m: _material_matches(m, material_select) if pd.notna(m) else False
            )
        ]

    df = df.dropna(subset=[x_col])
    if df.empty:
        print(f"No data with '{x_axis_var}' values found.")
        return plt.figure()

    # ---- Set up plot ----
    fig_size = set_plot_style(export_data=export_data)
    fig, ax = plt.subplots(figsize=fig_size)

    # Deduplicate temperatures that round to the same value
    # (different samples may have e.g. 301.8 and 302.0 → both round to 302)
    df['temperature'] = df['temperature'].round(0).astype(int)

    # Standard color cycle and marker list (same as rest of notebook)
    _cycle_colors = sns.color_palette("colorblind", 12)
    _cycle_markers = ['x', 'o', '*', 'd', '^', 'v', '+', '<', '>', 'p',
                      'P', 'h', 'H', 'X', 'D', '|', '_', '1', '2', '3',
                      '4', '8', 's']

    # Temperature color/marker assignment — uses the prop_cycle palette
    _dual_y = (y_axis_var_right is not None)
    unique_temps = sorted(df['temperature'].unique())
    temp_colors = {}
    temp_markers = {}
    for i, t in enumerate(unique_temps):
        if _dual_y:
            temp_colors[t] = 'black'  # dual y uses blue/red for variables
        else:
            temp_colors[t] = _cycle_colors[i % len(_cycle_colors)]
        temp_markers[t] = _cycle_markers[i % len(_cycle_markers)]

    # Y-axis label lookup (electrical + XRD)
    # Determine 2D vs 3D carrier density label from film_thickness
    _is_2d = (ppms_objects and ppms_objects[0].film_thickness == 1)
    _n_label_axis = r'$n \cdot t$ (cm$^{-2}$)' if _is_2d else r'$n$ (cm$^{-3}$)'
    _n_label_symbol = r'$n \cdot t$' if _is_2d else r'$n$'

    _y_labels = {
        # Electrical
        'mobility': r'$\mu$ (cm$^2$V$^{-1}$s$^{-1}$)',
        'carrier_density': _n_label_axis,
        'rho_xx': r'$\rho_{xx}$ ($\Omega \cdot$m)',
        'sheet_resistance': r'$R_{\square}$ ($\Omega/\square$)',
        'hall_coefficient': r'$R_H$ (m$^3$/C)',
        # XRD / structural
        'a': r'$a_{\parallel}$ $(\mathrm{\AA})$',
        'c': r'$c_{\perp}$ $(\mathrm{\AA})$',
        'a_and_c': r'Lattice parameter $(\mathrm{\AA})$',
        'qx': r'$Q_{\parallel}$ $(2\pi/\mathrm{\AA})$',
        'qz': r'$Q_{\perp}$ $(2\pi/\mathrm{\AA})$',
        'intensity': 'Fitted intensity',
        'intensity_measured': 'Measured intensity',
        'intensity_norm': 'Normalised intensity',
        'intensity_norm_perpulse': 'Normalised intensity per pulse',
        'sigma_x': r'$\sigma_x$',
        'sigma_z': r'$\sigma_z$',
        'eta': r'$\eta$ (mixing)',
        'background': 'Background',
        'tetragonality': r'$c/a$',
    }

    _x_labels = {
        'Ba_doping': r'Ba doping (\%)',
        'La_doping': r'La doping (\%)',
        'O2_growth_pressure': r'$\mathrm{O_2}$ growth pressure (mbar)',
        'T_growth': r'Growth temperature ($^{\circ}$C)',
        'n_pulses': r'Number of pulses',
        'thickness_nm': r'Thickness (nm)',
        'Hz': r'Repetition rate (Hz)',
        'E_mJ': r'Laser energy (mJ)',
    }

    # --- Helper to plot one variable on an axis ---
    def _plot_var(axis, var_name, df_in, use_temp_color=True, override_color=None):
        if var_name == 'a_and_c':
            if 'a' not in df_in.columns or 'c' not in df_in.columns:
                print("  Columns 'a' and/or 'c' not found in data — no XRD data available for this layer.")
                return
            sub_a = df_in.dropna(subset=[x_col, 'a']).sort_values(x_col)
            sub_c = df_in.dropna(subset=[x_col, 'c']).sort_values(x_col)
            if not sub_a.empty:
                axis.plot(sub_a[x_col].astype(float).values,
                          sub_a['a'].astype(float).values,
                          marker='o', linestyle='None', color='red', zorder=5)
            if not sub_c.empty:
                axis.plot(sub_c[x_col].astype(float).values,
                          sub_c['c'].astype(float).values,
                          marker='s', linestyle='None', color='blue', zorder=5)
        else:
            if var_name not in df_in.columns:
                print(f"  Column '{var_name}' not found in data.")
                return
            err_col = var_name + '_error'
            has_errors = err_col in df_in.columns
            for t in sorted(df_in['temperature'].unique()):
                df_t = df_in[df_in['temperature'] == t].dropna(subset=[x_col, var_name])
                if df_t.empty:
                    continue
                df_t = df_t.sort_values(x_col)
                x_vals = df_t[x_col].astype(float).values
                y_vals = df_t[var_name].astype(float).values
                mk = temp_markers.get(t, 'o')

                if override_color is not None:
                    c = override_color
                elif use_temp_color:
                    c = temp_colors.get(t, 'black')
                else:
                    c = 'black'

                if has_errors:
                    yerr = df_t[err_col].values
                    valid = pd.notna(yerr)
                    yerr_clean = np.where(valid, np.abs(yerr.astype(float)), 0)
                    if np.any(yerr_clean > 0):
                        axis.errorbar(x_vals, y_vals, yerr=yerr_clean,
                                      fmt=mk, color=c, linestyle='None', zorder=5)
                        continue
                axis.plot(x_vals, y_vals, marker=mk, linestyle='None',
                          color=c, zorder=5)

    # Determine left-axis colour
    if y_axis_var == 'a_and_c':
        left_color = None  # uses red/blue internally
    elif y_axis_var_right is not None:
        left_color = 'blue'
    else:
        left_color = None  # use temperature colormap

    use_temp_left = (left_color is None and y_axis_var != 'a_and_c')
    _plot_var(ax, y_axis_var, df, use_temp_color=use_temp_left,
              override_color=left_color)

    y_label = _y_labels.get(y_axis_var, y_axis_var.replace('_', ' ').title())
    ax.set_ylabel(y_label, color=left_color if left_color else 'black')
    if left_color:
        ax.tick_params(axis='y', labelcolor=left_color)

    if log_y:
        ax.set_yscale('log')
    if log_x:
        ax.set_xscale('log')

    # --- Right y-axis ---
    ax_right = None
    if y_axis_var_right is not None:
        right_color = 'red'
        ax_right = ax.twinx()
        _plot_var(ax_right, y_axis_var_right, df, use_temp_color=False,
                  override_color=right_color)
        ax_right.set_ylabel(
            _y_labels.get(y_axis_var_right, y_axis_var_right.replace('_', ' ').title()),
            color=right_color)
        # Disable grid on right y-axis to avoid double grid lines
        ax_right.grid(False)
        if y_lim_right is not None:
            ax_right.set_ylim(y_lim_right)
        if log_y_right:
            ax_right.set_yscale('log')
        # Apply tick colors AFTER log scale so new ticks inherit the color
        ax_right.tick_params(axis='y', labelcolor=right_color)

        # Add 10% padding to right y-axis so dual-axis traces don't overlap
        if y_lim_right is None:
            _lo, _hi = ax_right.get_ylim()
            if ax_right.get_yscale() == 'log':
                _ratio = _hi / _lo if _lo > 0 else 1
                _pad = _ratio ** 0.1
                ax_right.set_ylim(_lo / _pad, _hi * _pad)
            else:
                _span = _hi - _lo if _hi != _lo else abs(_hi) * 0.1 or 1
                ax_right.set_ylim(_lo - 0.1 * _span, _hi + 0.1 * _span)

    # --- Line plot and fit line ---
    def _line_and_fit(axis, var_name, color, df_in):
        if var_name == 'a_and_c':
            for v, c in [('a', 'red'), ('c', 'blue')]:
                sub = df_in.dropna(subset=[x_col, v])[[x_col, v]].copy()
                sub = sub.sort_values(x_col)
                xv = sub[x_col].astype(float).values
                yv = sub[v].astype(float).values
                if len(xv) < 2:
                    continue
                if line_plot:
                    axis.plot(xv, yv, color=c, zorder=4)
                if fit_line:
                    slope, intercept, r, _, _ = linregress(xv, yv)
                    x_fit = np.linspace(xv.min(), xv.max(), 200)
                    axis.plot(x_fit, slope * x_fit + intercept, color=c,
                              linestyle='--', zorder=4,
                              label=f'fit ($R^2$={r**2:.3f})')
        else:
            if var_name not in df_in.columns:
                return
            # For electrical vars, do line/fit per temperature
            for t in sorted(df_in['temperature'].unique()):
                df_t = df_in[df_in['temperature'] == t]
                sub = df_t.dropna(subset=[x_col, var_name])[[x_col, var_name]].copy()
                sub = sub.sort_values(x_col)
                xv = sub[x_col].astype(float).values
                yv = sub[var_name].astype(float).values
                if len(xv) < 2:
                    continue
                c = color if color else temp_colors.get(t, 'black')
                if line_plot:
                    axis.plot(xv, yv, color=c, zorder=4)
                if fit_line:
                    slope, intercept, r, _, _ = linregress(xv, yv)
                    x_fit = np.linspace(xv.min(), xv.max(), 200)
                    axis.plot(x_fit, slope * x_fit + intercept, color=c,
                              linestyle='--', zorder=4,
                              label=f'{t:.0f} K fit ($R^2$={r**2:.3f})')

    if line_plot or fit_line:
        _line_and_fit(ax, y_axis_var, left_color, df)
        if ax_right is not None and y_axis_var_right is not None:
            _line_and_fit(ax_right, y_axis_var_right, 'red', df)

    # Bulk reference lines for lattice parameters
    if show_bulk_refs and y_axis_var in ('a', 'c', 'a_and_c'):
        plotted_layer_mats = df['layer_material'].dropna().unique() if 'layer_material' in df.columns else []
        for key, mat_info in materials.items():
            if key == 'STO':
                continue
            if not any(_material_matches(key, [m]) or _material_matches(m, [key]) for m in plotted_layer_mats):
                continue
            ax.axhline(y=mat_info['a_bulk'], color='grey', linestyle='--',
                       linewidth=2.0, alpha=0.5, label=f'{key} bulk')

    # X-axis label
    if _input_new_mode:
        ax.set_xlabel(_input_new_label)
    elif _sample_name_mode:
        ax.set_xlabel('Sample')
        ax.set_xticks(range(len(_sample_order)))
        ax.set_xticklabels(_sample_order, rotation=45, ha='right', fontsize=7)
    else:
        ax.set_xlabel(_x_labels.get(x_axis_var, x_axis_var.replace('_', ' ').title()))

    if x_lim is not None:
        ax.set_xlim(x_lim)
    if y_lim is not None:
        ax.set_ylim(y_lim)

    # --- Legend ---
    legend_elements = []

    # Temperature entries — skip when using dual y-axis (colour encodes variable, not T)
    _xrd_only_vars = ('a', 'c', 'a_and_c', 'tetragonality',
                      'sigma_x', 'sigma_z', 'eta',
                      'background', 'intensity',
                      'intensity_norm', 'intensity_norm_perpulse')
    if not _dual_y:
        if len(unique_temps) > 1:
            for t in unique_temps:
                legend_elements.append(
                    Line2D([0], [0], marker=temp_markers[t], color=temp_colors[t],
                           markerfacecolor=temp_colors[t], markersize=4,
                           linestyle='None', label=f'{t:.0f} K'))
        elif len(unique_temps) == 1 and y_axis_var not in _xrd_only_vars:
            t0 = unique_temps[0]
            legend_elements.append(
                Line2D([0], [0], marker=temp_markers[t0], color=temp_colors[t0],
                       markerfacecolor=temp_colors[t0], markersize=4,
                       linestyle='None', label=f'{t0:.0f} K'))

    # Dual y-axis colour entries
    if y_axis_var == 'a_and_c':
        legend_elements.append(mpatches.Patch(facecolor='red', edgecolor='black',
                                              linewidth=0.5, label=r'$a_{\parallel}$ (in-plane)'))
        legend_elements.append(mpatches.Patch(facecolor='blue', edgecolor='black',
                                              linewidth=0.5, label=r'$c_{\perp}$ (out-of-plane)'))
    elif y_axis_var_right is not None:
        _legend_names = {
            'mobility': r'$\mu$',
            'carrier_density': _n_label_symbol,
            'rho_xx': r'$\rho_{xx}$',
            'sheet_resistance': r'$R_{\square}$',
            'hall_coefficient': r'$R_H$',
            'a': r'$a_{\parallel}$',
            'c': r'$c_{\perp}$',
            'tetragonality': r'$c/a$',
        }
        legend_elements.append(mpatches.Patch(facecolor='blue', edgecolor='black',
                                              linewidth=0.5, label=_legend_names.get(y_axis_var, y_axis_var.replace('_', ' ').title())))
        legend_elements.append(mpatches.Patch(facecolor='red', edgecolor='black',
                                              linewidth=0.5, label=_legend_names.get(y_axis_var_right, y_axis_var_right.replace('_', ' ').title())))

    if show_key is not False and legend_elements:
        loc = show_key if isinstance(show_key, str) else 'best'
        ax.legend(handles=legend_elements, loc=loc)

    fig.tight_layout()

    # Export
    if export_data:
        _x_label_export = (_input_new_label.replace(' ', '_') if _input_new_mode and _input_new_label
                           else x_axis_var)
        label = f'growth_opt_{_x_label_export}_vs_{y_axis_var}'
        if y_axis_var_right:
            label += f'_and_{y_axis_var_right}'
        # Try to save to the directory of the first PPMSData object
        first_ppms = _get_ppms_object(sample_list, sample_codes[0]) if sample_codes else None
        if first_ppms is not None and first_ppms.directory is not None:
            out_path = Path(first_ppms.directory) / f'{label}.{fig_format}'
            fig.savefig(str(out_path), dpi=600, bbox_inches='tight', transparent=True)
            print(f"  Saved: {out_path}")
        else:
            print("  export_data=True but no output directory found on first sample.")

    plt.show()
    return fig
