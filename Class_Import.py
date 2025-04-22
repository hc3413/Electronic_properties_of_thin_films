#import all the libraries needed
from import_dep import *

# Define a class to store the VNA data alongside its ascociated with the filename, device index, state and DC resistance
@dataclass
class PPMSData:
    data_import_np: np.ndarray = None #raw data from the PPMS instrument
    data_import_df: pd.DataFrame = None #data from the PPMS instrument in a pandas dataframe
    
    data_np: np.ndarray = None #PPMS data sliced and reordered for analysis
    data_np_nd: np.ndarray = None #PPMS data sliced into a multi-dimensional numpy array with the dimensions (temperature, field, current, columns, index)
    data_df: pd.DataFrame = None #PPMS data sliced and reordered for analysis in a pandas dataframe
    
    ctf: list = None # Storing extract values from data [current_unique, temp_unique, field_unique, current_no, temp_no, field_no]
    
    tf_av: np.ndarray = None #average temperature and field values as measured during each set of current measurements
    
    res_data: np.ndarray = None #resistivity data calculated using the Van der Pauw method or using a HallBar: columns = ['Temp (K)', 'Field (T)', 'rho_xx_A (ohm.m)', 'rho_xx_B(ohm.m)','rho_xx_average(ohm.m)']
    res_data_df: pd.DataFrame = None #resistivity data calculated using the Van der Pauw method or a HallBar in a pandas dataframe
    R_squared_res: list = None #R-squared values for the linear regression of the IV data
    
    mag_res: np.ndarray = None #magnetoresistance data calculated for each temperature and field strength
    
    hall_data: np.ndarray = None #Hall resistivity data calculated for each temperature and field strength: columns = ['Temp (K)', 'Field (T)', 'rho_xy_A(ohm.m)', 'R_squared(I)_A', 'rho_xy_B(ohm.m)','R_squared(I)_B', 'rho_xy_average(ohm.m)']
    hall_data_df: pd.DataFrame = None #Hall resistivity data calculated for each temperature and field strength in a pandas dataframe
    hall_coefficient: np.ndarray = None #Hall coefficient data calculated for each temperature: columns = ['Temp (K)', 'Hallco_A', 'R^2(H)_A', 'Hallco_B','R^2(H)_B', 'Hallco_average']
    hall_coefficient_df: pd.DataFrame = None #Hall coefficient data calculated for each temperature in a pandas dataframe
    
    directory: str = None #directory to save the output files to
    filename: str = None #filename of the data
    film_thickness: float = None  #thickness of the film in meters - set to 1 as default, if set to 1 then you are calculating sheet resistance not resistivity
    material: str = None #material of the film
    plot_str: str = None #strings to be used at the start of the figure file name, including the sample code
    sample_code: str = None # Sample code to identify the sample
    measurement_type: str = 'VDP' #measurement type to identify the measurement type as VDP (Van der Pauw Geometry) or HallBar (Hall Bar Geometry)

def import_ppms_data(
    path, #path to the directory containing the PPMS data
    film_thickness = 1.0, #thickness of the film in meters - set to 1 as default, if set to 1 then you are calculating sheet resistance not resistivity
    material = 'NA', #materials of the thin film stack
    sample_code = 'NA', #strings to be used at the start of the figure file name, including the sample code
    V_inv = False, #invert the sense voltage in case the cables were switched during the measurement (for HC011)
    measurement_type: str = 'VDP' #measurement type to identify the measurement type as VDP or HallBar
    ):
    '''
    This function imports the data from the PPMS instrument.
    Meausured in the Van der Pauw configuration with 6 indices - 0/1 are the Hall Voltage measurements 2/3/4/5 are the resistivity measurements.
    Meausred in the HallBar configuration with 4 indices - 0/1 are the Hall Voltage measurements 2/3 are the resistivity measurements.
    Note that it converts the field from Oe to Tesla. (strictly speaking from H to B)
    '''
    PPMS_files = []
    # Get a list of all the files in the directory converting them to Path objects
    files = [Path(f) for f in os.listdir(path) if not f.startswith('.')]
    # Sort the files alphabetically
    files.sort()
    # Loop over all files in the directory extracting the data as a df, converting it to a numpy array, and slicing it to separete index values into a third dimension
    for count, fi in enumerate(files):  
                    
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
        elif measurement_type == 'HallBar':
            slice_size = 4
        else:
            raise ValueError("Invalid measurement type. Use 'VDP' or 'HallBar'.")

        # Slice the array, then stack these 2d arrays to generate a 3d array
        ppms_np = np.array([ndnp[q::6,:-1] for q in range(6)]) # the :-1 is to remove the index column
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
        PPMS_files.append(PPMSData(data_import_np = ppms_np, data_import_df = ppms_df, filename = fi, film_thickness = film_thickness, material = material, directory = path_out, plot_str = fi, sample_code = sample_code))
        
        # Print the filenames of the imported data to check they are correct along with the shape of the numpy array
        print(f'File {count+1} imported: {fi} with shape {ppms_np.shape}')
    return PPMS_files



def import_all_datasets():
    """
    Import all electrical measurement datasets from various samples.
    Separately imports the VDP and HallBar datasets and returns a list of PPMSData objects for each
    
    Returns:
        list: A list of PPMSData objects containing all imported datasets.
    """
    all_directories = []
    dat_raw = []
    
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
    ]
    
    # Import each dataset
    for i, dataset in enumerate(datasets):
        path = Path(dataset['path'])
        all_directories.append(path)
        
        # Extract V_inv parameter if it exists, otherwise default to False
        v_inv = dataset.get('V_inv', False)
        # Extract the measurement type if it exists, otherwise default to 'VDP'
        measurement_type = dataset.get('measurement_type', 'VDP')
        
        # Import the data
        dat_raw.append(import_ppms_data(
            path,
            film_thickness=dataset['film_thickness'],
            material=dataset['material'],
            sample_code=dataset['sample_code'],
            V_inv=v_inv,
            measurement_type=measurement_type
        ))
    
    print(f"Successfully imported {len(datasets)} datasets.")
    return dat_raw
