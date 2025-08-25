#import all the libraries needed
from import_dep import *

# Define a class to store the VNA data alongside its ascociated with the filename, device index, state and DC resistance
@dataclass
class PPMSData:
    # Raw Data
    data_import_np: np.ndarray = None #raw data from the PPMS instrument
    data_import_df: pd.DataFrame = None #data from the PPMS instrument in a pandas dataframe
    
    # Sorted Data
    data_np: np.ndarray = None #PPMS data sliced and reordered for analysis
    data_np_nd: np.ndarray = None #PPMS data sliced into a multi-dimensional numpy array with the dimensions (temperature, field, current, columns, index)
    data_df: pd.DataFrame = None #PPMS data sliced and reordered for analysis in a pandas dataframe
    
    # Storage Parameters for field/current/temperature values
    ctf: list = None # Storing extract values from data [current_unique, temp_unique, field_unique, current_no, temp_no, field_no]
    tf_av: np.ndarray = None #average temperature and field values as measured during each set of current measurements
    
    # Calculated Resistivity data
    res_data: np.ndarray = None # Bulk resistivity data (rho_xx) calculated using VDP or HallBar: columns = ['Temp (K)', 'Field (T)', 'rho_xx_A (Ohm.m)', 'rho_xx_B (Ohm.m)','rho_xx_average (Ohm.m)', 'rho_error (Ohm.m)']
    res_data_df: pd.DataFrame = None # Bulk resistivity data (rho_xx) calculated using VDP or HallBar in a pandas dataframe
    R_squared_res: list = None #R-squared values for the linear regression of the I-V data used for rho_xx calculation
    mag_res: np.ndarray = None #magnetoresistance data calculated for each temperature and field strength: columns = [MR_A(%), MR_B(%), MR_avg(%), MR_avg_error(%)]
    
    # Calculated Hall data
    hall_data: np.ndarray = None #Hall resistivity data (rho_xy) calculated for each temperature and field strength: columns = ['Temp (K)', 'Field (T)', 'rho_xy_A (Ohm.m)', 'R^2(I-V)_A', 'rho_xy_B (Ohm.m)','R^2(I-V)_B', 'rho_xy_average (Ohm.m)']
    hall_data_df: pd.DataFrame = None #Hall resistivity data (rho_xy) calculated for each temperature and field strength in a pandas dataframe
    hall_coefficient: np.ndarray = None #Hall coefficient (R_H) and derived data calculated for each temperature: columns = ['Temp (K)', 'HallCoeff_A (m3/C)', 'R^2(B)_A', 'HallCoeff_B (m3/C)','R^2(B)_B', 'HallCoeff_Avg (m3/C)','R^2(B)_Avg', 'n (cm^-3)', 'n_Error (cm^-3)', 'Mobility (cm^2/Vs)', 'Mobility_Error (cm^2/Vs)']
    hall_coefficient_df: pd.DataFrame = None #Hall coefficient (R_H) and derived data calculated for each temperature in a pandas dataframe
    
    # Fitted Data



    # File and directory parameters
    directory: str = None #directory to save the output files to
    filename: str = None #filename of the data
    film_thickness: float = None  #thickness of the film in meters - set to 1 as default, if set to 1 then you are calculating sheet resistance not resistivity
    material: str = None #material of the film
    plot_str: str = None #strings to be used at the start of the figure file name, including the sample code
    sample_code: str = None # Sample code to identify the sample
    
    # Measurement parameters
    measurement_type: str = 'VDP' #measurement type to identify the measurement type as VDP (Van der Pauw Geometry) or HallBar (Hall Bar Geometry)
    hb_dimensions: tuple = None # dimensions of the Hall Bar in um: (width of channel, length between voltage arms)
    rotator: bool = False #If True, rotator used, theta angle is stored in the V_source column and used to calculate the field strength

def import_ppms_data(
    path, #path to the directory containing the PPMS data
    film_thickness = 1.0, #thickness of the film in meters - set to 1 as default, if set to 1 then you are calculating sheet resistance not resistivity
    material = 'NA', #materials of the thin film stack
    sample_code = 'NA', #strings to be used at the start of the figure file name, including the sample code
    V_inv: bool = False, #invert the sense voltage in case the cables were switched during the measurement (for HC011)
    measurement_type: str = 'VDP', #measurement type to identify the measurement type as 'VDP' or 'HallBar'
    hb_dimensions: tuple = None, # Hall bar dimensions in um (width, length_between_arms), required if measurement_type='HallBar'
    
    ):
    '''
    This function imports the data from the PPMS instrument.
    Measured in the Van der Pauw configuration with 6 indices - 0/1 are the Hall Voltage measurements 2/3/4/5 are the resistivity measurements.
    Measured in the HallBar configuration with 4 indices - 0/1 are the resistivity measurements 2/3 are the Hall Voltage measurements.
    Note that it converts the field from Oe to Tesla. (strictly speaking from H to B)
    '''
    PPMS_files = []
    # Get a list of all the files in the directory converting them to Path objects
    files = [Path(f) for f in os.listdir(path) if not f.startswith('.')]
    # Sort the files alphabetically
    files.sort()
    # Loop over all files in the directory extracting the data as a df, converting it to a numpy array, and slicing it to separete index values into a third dimension
    for count, fi in enumerate(files):  
        
        
        # --- Determine if rotator was used by checking headers ---
        rotator = False # Default to False
        try:
            # Read the header row (row 5, which is skiprows=4 since it's 0-indexed)
            header_df = pd.read_csv(path.joinpath(fi), sep='\t', skiprows=4, nrows=1, comment='N')
            # Check if 'Angle (deg)' is present in any column header
            if any('Angle (deg)' in col for col in header_df.columns):
                rotator = True
                print(f"Rotator detected in file: {fi.name}")
        except Exception as e:
            print(f"Could not read headers to check for rotator in file: {fi}, {e}")
        # --- End rotator check ---
        
                    
        # Import the file, treating lines with non-numeric data as comments
        try:
            ppms_df = pd.read_csv(path.joinpath(fi), sep='\t', skiprows=6, header=None, comment='N')
        except Exception as e:
            print(f"Error with file: {fi}, {e}")
            continue

        # Drop unwanted columns
        try:
            ppms_df.drop(ppms_df.columns[[6, 7]], axis=1, inplace=True)
        except Exception as e:
            print(f"Error with file: {fi}, {e}")
            continue
        
    
        

        # Assign headers to columns
        ppms_df.columns = ['Temp (K)', 'Field', 'Source (A)', 'Source (V)', 'Sense (V)', 'index']

        # Convert columns to the appropriate data types
        ppms_df['Temp (K)'] = pd.to_numeric(ppms_df['Temp (K)'], errors='coerce')
        ppms_df['Field'] = pd.to_numeric(ppms_df['Field'], errors='coerce')
        ppms_df['Source (A)'] = pd.to_numeric(ppms_df['Source (A)'], errors='coerce')
        ppms_df['Source (V)'] = pd.to_numeric(ppms_df['Source (V)'], errors='coerce')
        ppms_df['Sense (V)'] = pd.to_numeric(ppms_df['Sense (V)'], errors='coerce')
        ppms_df['index'] = pd.to_numeric(ppms_df['index'], downcast='integer', errors='coerce')
        
        # If the sample holder is rotated in the field, the theta angle is stored in the 'Source (V)' column
        # The field component perpendicular to the sample surface is then calculated using the theta angle
        if rotator == True:
            # Rename the 'Source (V)' column to 'Theta (degrees)' and re-write the field column using the theta angle to get the perpendicular field
            ppms_df.rename(columns={'Source (V)': 'Theta (degrees)'}, inplace=True)
            # Convert the theta angle from degrees to radians
            ppms_df['Theta (degrees)'] = np.deg2rad(ppms_df['Theta (degrees)'])
            # Calculate the field using the theta angle
            ppms_df['Field'] = (ppms_df['Field'] * np.cos(ppms_df['Theta (degrees)']))
    
        if V_inv == True:
            # Multiply the 'Sense (V)' column by -1 for only HC011 as the cables were switched the wrong way around
            ppms_df['Sense (V)'] *= -1
        
        # Drop rows where all values are NaN (e.g., empty lines)
        ppms_df.dropna(how='all', inplace=True)

        # Convert the 'index' column to integer format
        ppms_df['index'] = ppms_df['index'].astype(int)
        
        # Convert the field from Oe to Tesla
        ppms_df['Field'] = ppms_df['Field'] / 10000

        # Change the header to 'Field (T)'
        ppms_df.rename(columns={'Field': 'Field (T)'}, inplace=True)

        # Convert the data frame to a numpy array
        ndnp = ppms_df.to_numpy()
        
        # Slice the array taking every 6th or 4th value (i.e the different index values) 
        # 6th for the VDP data and 4th for the HallBar data
        if measurement_type == 'VDP':
            slice_size = 6
            indexes = 6
        elif measurement_type == 'HallBar':
            slice_size = 4
            indexes = 4
            if hb_dimensions is None:
                 raise ValueError("hb_dimensions must be provided for measurement_type='HallBar'")
        else:
            raise ValueError("Invalid measurement type. Use 'VDP' or 'HallBar'.")

        # Slice the array, then stack these 2d arrays to generate a 3d array
        # The :-1 is to remove the original index column from the raw data
        ppms_np = np.array([ndnp[q::slice_size,:-1] for q in range(slice_size)]) 
        # Re-order the array to have the shape (rows, columns, index)
        ppms_np = ppms_np.transpose(1,2,0)   

        
        ### Define the output directory to save plots and PowerPoint to - always defined by the first data set
        # Ensure the path is a string and replace '/Data/' with '/Output/'
        path_str = str(path)
        path_out_str = path_str.replace('/Data', '/Output')
        
        # Ensure the path ends with a trailing slash
        if not path_out_str.endswith('/'):
            path_out = Path(path_out_str + '/')
            
        # Convert back to Path object
        path_out = Path(path_out_str)
        
        # Store the extracted data into a list of PPMSData objects
        PPMS_files.append(PPMSData(
            data_import_np = ppms_np, 
            data_import_df = ppms_df, 
            filename = fi.name, # Use fi.name to get just the filename string
            film_thickness = film_thickness, 
            material = material, 
            directory = path_out, 
            plot_str = fi.stem, # Use fi.stem for filename without extension
            sample_code = sample_code,
            measurement_type = measurement_type,
            hb_dimensions = hb_dimensions, # Store dimensions if HallBar
            rotator = rotator # Whether the rotator was used
            ))
        
        # Print the filenames of the imported data to check they are correct along with the shape of the numpy array
        print(f'File {count+1} imported: {fi.name} with shape {ppms_np.shape}')
    return PPMS_files



def import_all_datasets():
    """
    Import all electrical measurement datasets from various samples.
    Separately imports the VDP and HallBar datasets and returns a list of PPMSData objects for each
    
    Returns:
        list: A list of PPMSData objects containing all imported datasets.
        all_directories: A list of directories where the datasets are stored which is used for both VDP and HallBar datasets but the type is specified
        There are separate output lists for VDP and HallBar datasets so they can have different processing functions applied to them
        
        'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_SR002_M_WO3_SiO2_Si/Data/',
        'film_thickness': 20e-9,
        'material': '$WO_3/SiO_2/Si$',
        'sample_code': 'SR002_M',
        'V_inv': False,
        'notes': 'medium hydrogen concentration annealed sample'
        'measurement_type' = ('VDP' or 'HallBar'): A string indicating the type of measurement 
        'hb_dimensions' = (width of channel, length between voltage arms): dimensions of the Hall Bar in um - set to None by default
        'rotator' = True or False: A boolean indicating if the rotator was used in the measurement
                    If rotator is True then the theta angle is stored in the V_source column and used to calculate the field strength
                    This can be done as the rotator is only ever used with the K6221 which doesn't have a source voltage output
       
    """
    all_directories = []
    dat_raw_vdp = []
    dat_raw_hb = []
    
    # Define all datasets as a list of dictionaries for easier maintenance
    datasets = [
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_JT029_SrRuO3/Data/',
            'film_thickness': 14e-9,
            'material': '$SrRuO_3/SrTiO_3$',
            'sample_code': 'JT029',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_HC003_LaScO3/Data/',
            'film_thickness': 1,
            'material': '$LaScO_3/SrTiO_3$',
            'sample_code': 'HC003',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_YFW042_HZO/Data/',
            'film_thickness': 1,
            'material': '$HZO$',
            'sample_code': 'YFW042',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_DR004_BSO_LSO/Data/',
            'film_thickness': 1,
            'material': '$BaSnO_3/LaScO_3/SrTiO_3$',
            'sample_code': 'DR004',
            'notes': 'BSO layer is 309 nm, LSO layer is 12.9 nm'
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_DR003_BSO_LSO/Data/',
            'film_thickness': 1,
            'material': '$BaSnO_3/LaScO_3/SrTiO_3$',
            'sample_code': 'DR003',
            'notes': 'BSO layer is 317 nm, LSO layer is 12.7 nm'
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_HC006_SSO/Data/',
            'film_thickness': 1,
            'material': '$La:SrSnO_3/SrTiO_3$',
            'sample_code': 'HC006',
            'notes': 'BSO layer is 317 nm, LSO layer is 12.7 nm'
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_DR001_BSO_LSO/Data/',
            'film_thickness': 1,
            'material': '$BaSnO_3/LaScO_3/SrTiO_3$',
            'sample_code': 'DR001',
            'notes': 'BSO layer is 317 nm, LSO layer is 12.7 nm'
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_YFW045_HZO/Data/',
            'film_thickness': 1,
            'material': '$HZO$',
            'sample_code': 'YFW045',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_YFW044_HZO/Data/',
            'film_thickness': 1,
            'material': '$HZO$',
            'sample_code': 'YFW044',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_YFW046_HZO/Data/',
            'film_thickness': 1,
            'material': '$HZO$',
            'sample_code': 'YFW046',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_HC011_BSO_LSO_LaSSO_STO/Data/',
            'film_thickness': 1,
            'material': '$BaSnO_3/LaScO_3/La:SrSnO_3/SrTiO_3$',
            'sample_code': 'HC011',
            'V_inv': True,
            'notes': 'Using V_inv = True to invert the voltage data for this set as the cables were reversed'
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_HC013_LSO_BSO_LaSSO_STO/Data/',
            'film_thickness': 1,
            'material': '$LaScO_3/BaSnO_3/La:SrSnO_3/SrTiO_3$',
            'sample_code': 'HC013',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_HC014_BSO_LSO_LaSSO_STO/Data/',
            'film_thickness': 1,
            'material': '$BaSnO_3/LaScO_3/La:SrSnO_3/SrTiO_3$',
            'sample_code': 'HC014',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_DV127_WO3_YAO/Data/',
            'film_thickness': 140e-9,
            'material': '$WO_3/YAlO_3$',
            'sample_code': 'DV127',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_DR008a_BSO_LSO_STO/Data/',
            'film_thickness': 1,
            'material': '$LaScO_3/BaSnO_3/SrTiO_3$',
            'sample_code': 'DR008a',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_HC016_LSO_BSO_LaSSO_STO/Data/',
            'film_thickness': 1,
            'material': '$LaScO_3/BaSnO_3/La:SrSnO_3/SrTiO_3$',
            'sample_code': 'HC016',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_HC015_LSO_BSO_STO/Data/',
            'film_thickness': 1,
            'material': '$LaScO_3/BaSnO_3/SrTiO_3$',
            'sample_code': 'HC015',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_DR009_BSO_LSO_STO/Data/',
            'film_thickness': 1,
            'material': '$LaScO_3/BaSnO_3/SrTiO_3$',
            'sample_code': 'DR009',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_DR010_BSO_LSO_STO/Data/',
            'film_thickness': 1,
            'material': '$LaScO_3/BaSnO_3/SrTiO_3$',
            'sample_code': 'DR010',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_DR012_LSO_BSO_BSO_STO/Data/',
            'film_thickness': 1,
            'material': '$LaScO_3/BaSnO_3/BaSnO_3/SrTiO_3$',
            'sample_code': 'DR012',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_JT082_BSO_BTO_SSTO_STO/Data/',
            'film_thickness': 1,
            'material': '$BaSnO_3/BaTiO_3/La_ySn_xSr_{1-x}TiO_3/SrTiO_3$',
            'sample_code': 'JT082',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_DR011_LSO_BSO_STO/Data/',
            'film_thickness': 1,
            'material': '$LaScO_3/BaSnO_3/SrTiO_3$',
            'sample_code': 'DR011',
            'notes': ''
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_HC011_BSO_LSO_LaSSO_STO/Data/',
            'film_thickness': 1,
            'material': '$BaSnO_3/LaScO_3/La:SrSnO_3/SrTiO_3$',
            'sample_code': 'HC011',
            'V_inv': False,
            'notes': 'duplicated and now not using Vin for later data'
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_SR003_S_WO3_SiO2_Si/Data/',
            'film_thickness': 20e-9,
            'material': '$WO_3/SiO_2/Si$',
            'sample_code': 'SR003_S',
            'V_inv': False,
            'notes': 'highest hydrogen concentration annealed sample'
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_SR002_M_WO3_SiO2_Si/Data/',
            'film_thickness': 20e-9,
            'material': '$WO_3/SiO_2/Si$',
            'sample_code': 'SR002_M',
            'V_inv': False,
            'notes': 'medium hydrogen concentration annealed sample'
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/HallBar1_DV127_WO3_YAO/Data/',
            'film_thickness': 140e-9,
            'material': '$WO_3/YAlO_3/STO$',
            'sample_code': 'DV127HB1',
            'V_inv': False,
            'measurement_type': 'HallBar',
            'hb_dimensions': (30, 100),  # width of channel, length between voltage arms in um
            'notes': 'processed into hall bar to look at disparity between different directions'
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/HallBar4_DV127_WO3_YAO/Data/',
            'film_thickness': 140e-9,
            'material': '$WO_3/YAlO_3/STO$',
            'sample_code': 'DV127HB4',
            'V_inv': False,
            'measurement_type': 'HallBar',
            'hb_dimensions': (30, 100),  # width of channel, length between voltage arms in um
            'notes': 'processed into hall bar to look at disparity between different directions - this sample has some damage to the electrodes where the wirebonding was done'
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_HC017_BSO_LSO_SSO_STO/Data/',
            'film_thickness': 1,
            'material': '$BaSnO_3/LaScO_3/La:SrSnO_3/SrTiO_3$',
            'sample_code': 'HC017',
            'V_inv': False,
            'notes': 'low oxygen growth pressure 0.01 for BSO'
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_HC018_BSO_LSO_SSO_STO/Data/',
            'film_thickness': 1,
            'material': '$BaSnO_3/LaScO_3/La:SrSnO_3/SrTiO_3$',
            'sample_code': 'HC018',
            'V_inv': False,
            'notes': 'high oxygen growth pressure 0.13 for BSO'
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_HP05_MoS2/Data/',
            'film_thickness': 2.5e-3,
            'material': '$MoS_2$',
            'sample_code': 'HP05',
            'V_inv': False,
            'notes': 'first Xin sample of MoS2 testing with silver paste contacts which were wirebonded on to and using quartz to mount it thermally'
        },
        {
            'path': '/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_YL196_ZO_STO/Data/',
            'film_thickness': 1,
            'material': '$ZrO_2/STO$',
            'sample_code': 'YL196',
            'V_inv': False,
            'notes': '2DEG at STO interface at low T?'
        },
    ]
    
    # Import each dataset
    for i, dataset in enumerate(datasets):
        path = Path(dataset['path'])
        all_directories.append(path)
        
        # Extract V_inv parameter if it exists, otherwise default to False
        v_inv = dataset.get('V_inv', False)
        
        # Extract the measurement type if it exists, otherwise default to 'VDP'
        measurement_type = dataset.get('measurement_type', 'VDP')
        
        # Extract Hall bar dimensions if measurement_type is 'HallBar'
        hb_dimensions = dataset.get('hb_dimensions', None) if measurement_type == 'HallBar' else None
        
  
    

        
        if measurement_type == 'VDP':
            # Import the VDP data
            dat_raw_vdp.append(import_ppms_data(
                path,
                film_thickness=dataset['film_thickness'],
                material=dataset['material'],
                sample_code=dataset['sample_code'],
                V_inv=v_inv,
                measurement_type=measurement_type,
                # No hb_dimensions needed for VDP

            ))
            
        elif measurement_type == 'HallBar':
            # Check if hb_dimensions were provided for HallBar type
            if hb_dimensions is None:
                print(f"Warning: hb_dimensions not provided for HallBar sample {dataset['sample_code']}. Skipping this dataset.")
                continue 

            # Import the HallBar data, passing hb_dimensions
            dat_raw_hb.append(import_ppms_data(
                path,
                film_thickness=dataset['film_thickness'],
                material=dataset['material'],
                sample_code=dataset['sample_code'],
                V_inv=v_inv,
                measurement_type=measurement_type,
                hb_dimensions=hb_dimensions # Pass the dimensions here

            ))

    
    print(f"Successfully imported {len(dat_raw_vdp)} VDP dataset(s) and {len(dat_raw_hb)} HallBar dataset(s).")
    return dat_raw_vdp, dat_raw_hb


