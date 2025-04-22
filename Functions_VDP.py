#import all the libraries needed
from import_dep import * # Import Dependencies
from Class_Import import * # Import the data storage class
from Functions_General import * # Import the general functions used to process the data

### ----------- Functions for extraction and analysis of VDP data ----------- ###

def vdp_equation(rho_sheet, R_A, R_B):
    '''Solves the Van der Pauw equation for the sheet resistivity of a thin film sample.'''
    return np.exp(-np.pi * R_A / rho_sheet) + np.exp(-np.pi * R_B / rho_sheet) - 1

def Error_Rs(R_sheet,R_32_10, R_20_31):
    '''Calculates the error in each calculation of sheet resistivity with the VDP formula by using error propogation
    R_sheet_A is the sheet resistivity calculated from the VDP formula
    f(R_sheet;Ra,Rb) = e^(-pi*Ra/R_sheet) + e^(-pi*Rb/R_sheet) - 1 = 0
    df/dRa = -pi/R_sheet * e^(-pi*Ra/R_sheet)
    df/dRb = -pi/R_sheet * e^(-pi*Rb/R_sheet)
    df/dR_sheet = pi/R_sheet^2 * (Ra*e^(-pi*Ra/R_sheet) + Rb*e^(-pi*Rb/R_sheet))
    dRsdRa = dR_sheet/dRa = -df/dRa/df/dR_sheet = R_sheet*e^(-pi*Ra/R_sheet)/(Ra*e^(-pi*Ra/R_sheet) + Rb*e^(-pi*Rb/R_sheet))
    dRsdRb = dR_sheet/dRb = -df/dRb/df/dR_sheet = R_sheet*e^(-pi*Rb/R_sheet)/(Ra*e^(-pi*Ra/R_sheet) + Rb*e^(-pi*Rb/R_sheet))'''
    R_A = R_32_10[0]
    dRa = R_32_10[4] #R_A_error
    R_B = R_20_31[0]
    dRb = R_20_31[4] #R_B_error
    # Error propogation for the sheet resistivity
    dRsdRa = (np.exp(-np.pi * R_A / R_sheet) * R_sheet)/(R_A*np.exp(-np.pi * R_A / R_sheet) + R_B*np.exp(-np.pi * R_B / R_sheet))
    dRsdRb = (np.exp(-np.pi * R_B / R_sheet) * R_sheet)/(R_A*np.exp(-np.pi * R_A / R_sheet) + R_B*np.exp(-np.pi * R_B / R_sheet))
    
    R_sheet_error = np.sqrt((dRsdRa*dRa)**2 + (dRsdRb*dRb)**2)
    
    return R_sheet_error
    
    


def vdp_resistivity(
    PPMS_files, #list of PPMSData objects
    print_val = False, # print the initial guess vs calculated resistivity
    resistivity_guess = 0, #initial guess for the sheet resistivity, if set to 0 then the initial guess is calculated from an approximation with scaling factor
    filt_kern = 0, # median filter kernel size, if set to 0 then no filter is applied
    filt_sigma = 0, # gaussian filter sigma, if set to 0 then no filter is applied
    threshold = 0 # z-score threshold for outlier detection, if set to 0 then no filter is applied, (lower threshold means more points are considered outliers (2 is typical value))
    ):
    '''Calculates the resistivity of a thin film sample using the Van der Pauw method.
    The measurments for this were done in index values of 2,3,4,5 
    Each index value corresponds to a different configuration of the source and sense leads.
    These 4 index configurations are split into 2 pairs of configurations, with each pair being able to generate the resistivity of the film independantly.
    Four configurations were used for robustness, to enable a comparison of the resistivity values obtained with the source/sense positions swapped.
    '''
    
    for ppms in PPMS_files:
        # Extract the required data from the PPMSData object
        film_thickness = ppms.film_thickness
        ctf = ppms.ctf
        tf_av = ppms.tf_av
        data_np = ppms.data_np
        data_np_nd = ppms.data_np_nd
        
        
        # Initialise a list to store the R^2 value for the IV data and check we have a good fit for the linear regression
        R_squared = []
        
        
        # Initialize an empty np aray with indices: (temp_index, field_index, data_colums) 
        # storing each temperature, field, and the corresponding resitivities for: config A, config B, average of A and B along with the error
        res_data = np.zeros((ctf[4],ctf[5], 6))

        #Loop over each temperature and field combination, calculating the sheet resistivity using the Van der Pauw method
        for Ti in range(ctf[4]): #for each temperature index
            for Bi in range(ctf[5]): #for each field index
                
                ##### Step 1: Using linear regression on the current-voltage data
                #R_ij_kl[0] = the slope(resistance), R_ij_kl[1] = intercept, R_ij_kl[2] = R-squared value

                #First pair of Van der Pauw configurations
                R_32_10 = linregress(data_np_nd[Ti,Bi,:,2,2], data_np_nd[Ti,Bi,:,4,2])    
                R_20_31 = linregress(data_np_nd[Ti,Bi,:,2,3], data_np_nd[Ti,Bi,:,4,3])
                #Second pair of Van der Pauw configurations
                R_01_23 = linregress(data_np_nd[Ti,Bi,:,2,4], data_np_nd[Ti,Bi,:,4,4])
                R_13_02 = linregress(data_np_nd[Ti,Bi,:,2,5], data_np_nd[Ti,Bi,:,4,5])

                # Append the R-squared value to the list
                R_squared.extend([R_32_10[2], R_20_31[2], R_01_23[2], R_13_02[2]])
                
                ### Initial guess for rho_sheet based off an approximate scaling factor from source I, sense V data
                if resistivity_guess == 0:
                    R_to_rho_scaling =  4 #scaling factor 
                    initial_guess = R_32_10[0]* R_to_rho_scaling 
                else:
                    # Use the user defined initial guess
                    initial_guess = resistivity_guess

                ##### Step 2: solve the Van der Pauw equation for both pairs of configurations

                # Solve for rho_sheet for the first pair (R_A)
                R_sheet_A = fsolve(vdp_equation, initial_guess, args=(R_32_10[0], R_20_31[0]))[0]

                # Solve for rho_sheet for the second pair (R_B)
                R_sheet_B = fsolve(vdp_equation, initial_guess, args=(R_01_23[0], R_13_02[0]))[0]
                
                # Solve for isotropic film using parallel R_A
                #R_sheet_A = fsolve(vdp_equation, initial_guess, args=(R_32_10[0], R_01_23[0]))[0]

                # Solve for isotropic film using parallel R_B
                #R_sheet_B = fsolve(vdp_equation, initial_guess, args=(R_20_31[0], R_13_02[0]))[0]
                
                #### Step 3: Calculate the error in the sheet resistivity for each pair of configurations
                # Calculate the error in the sheet resistivity for each pair of configurations
                R_sheet_A_error = Error_Rs(R_sheet_A,R_32_10, R_20_31)
                R_sheet_B_error = Error_Rs(R_sheet_B,R_01_23, R_13_02)
                

                #### Step 4: Calculate the average sheet resistivity and the error in the average
                
                # Average the two solutions for the final sheet resistivity
                R_sheet = (R_sheet_A + R_sheet_B) / 2
                
                # Option to print the initial guess and the calculated sheet resistivity to check the values are close
                if print_val == True:
                    print(f'(R_s_Guess, R_s_calc) = ({initial_guess:.1e}, {R_sheet:.1e}) ohm/sq - {ppms.filename}')
        
                # Calculate the error in the average sheet resistivity by progating the errors in the individual sheet resistivity values
                R_sheet_error = 0.5*np.sqrt(R_sheet_A_error**2 + R_sheet_B_error**2)
                resistivity_error = R_sheet_error*film_thickness # if the film thickness is 1 then this is the error in the sheet resistance
                
                #### Step 5: Insert the new row to the np data array using tf_av for temperature and field values that are averaged over all the currents
                res_data[Ti,Bi,:] = [tf_av[Ti,Bi,0], tf_av[Ti,Bi,1], R_sheet_A*film_thickness, R_sheet_B*film_thickness,R_sheet*film_thickness, resistivity_error]

            # Option to filter the resistivity vs field data
            if filt_kern != 0 or filt_sigma != 0 or threshold != 0:
                res_data[Ti,:,4] = filter_data(res_data[Ti,:,4], filt_kern, filt_sigma, threshold)
            
            
        # Flatten the res_data array to a 2D array so it can be put into a df for debugging
        res_data_flat = res_data.reshape((ctf[4]*ctf[5],6))  
        # Convert the numpy array to a pandas dataframe 
        res_data_df = pd.DataFrame(res_data_flat, columns=['Temp (K)', 'Field (T)', 'rho_xx_A (ohm.m)', 'rho_xx_B(ohm.m)','rho_xx_average(ohm.m)', 'rho_error(ohm.m)'])
        
        # Store the data in the PPMSData object
        ppms.res_data = res_data
        ppms.res_data_df = res_data_df
        ppms.R_squared_res = R_squared
        
        # print(f'res_data({ppms.filename})')
        # print(res_data_df.head())
        
    return PPMS_files






def vdp_hall(
    PPMS_files, # list of PPMSData objects
    filt_kern = 0, # median filter kernel size, if set to 0 then no filter is applied
    filt_sigma = 0, # gaussian filter sigma, if set to 0 then no filter is applied
    threshold = 0 # z-score threshold for outlier detection, if set to 0 then no filter is applied, (lower threshold means more points are considered outliers (2 is typical value))             
    ):
    '''Calculates the Hall resistivity at every temperature and field strength
    Uses linear regression on the current-voltage data to obtain the Hall resistivity
    The Hall resistivity is calculated for each configuration of the source and sense leads and the two are averaged
    Though caution is needed on this average as the Hall resistivity is a vector quantity and the two configurations may not be identical
    (rho_xy = thickness*V(hall)/I(source))'''
    #ctf = [current_unique, temp_unique, field_unique, current_no, temp_no, field_no]
    
    for ppms in PPMS_files:
        # Extract the required data from the PPMSData object
        film_thickness = ppms.film_thickness
        ctf = ppms.ctf
        tf_av = ppms.tf_av
        data_np_nd = ppms.data_np_nd
        
        # Initialize an empty np aray to store the temperature, field, hall resistance: in config A, R^2 for that fit, Hall resistance in  config B, R^2 for that, and Average 
        hall_data = np.zeros((ctf[4],ctf[5], 7))
        
        # Initialize an empty np array to store the Temperature, Hall coefficient A, R^2 A, Hall Coefficient B, and average Hall coefficients
        hall_coefficient = np.zeros((ctf[4], 11))

        #Loop over each temperature using regression on the hall_resistivity-field  data to obtain the Hall coefficient at each temperature
        for Ti in range(ctf[4]):
            #Loop over each field using regression on the current-voltage data to obtain the Hall resistivity at each field
            for Bi in range(ctf[5]):
                
                #### Step 1: Use linear regression on the current-voltage data to obtain the Hall resistivity at each field
                # R_ij_kl[0] = the slope(resistance), R_ij_kl[1] = intercept, R_ij_kl[2] = R-squared value     
                
                # Hall Resitance in configuration A
                R_13_42 = linregress(data_np_nd[Ti,Bi,:,2,0], data_np_nd[Ti,Bi,:,4,0]) #Regression on (x,y) = (I_source, V_measured) 
                # Hall resistance in configuration B (with source and sense leads swapped)
                R_24_31 = linregress(data_np_nd[Ti,Bi,:,2,1], data_np_nd[Ti,Bi,:,4,1]) #Regression on (x,y) = (H_applied, V_measured)
                
                # Average the two solutions for the final sheet resistivity
                R_hall_average = (R_13_42[0] + R_24_31[0]) / 2
                
                # Insert the new row to the np data array
                hall_data[Ti,Bi,:] = [tf_av[Ti,Bi,0], tf_av[Ti,Bi,1], R_13_42[0]*film_thickness, R_13_42[2], R_24_31[0]*film_thickness, R_24_31[2], R_hall_average*film_thickness]

            # Option to filter Data Hall resistivity vs field
            if filt_kern != 0 or filt_sigma != 0 or threshold != 0:
                hall_data[Ti,:,6] = filter_data(hall_data[Ti,:,6], filt_kern, filt_sigma, threshold)
            
                
            # Step 2: Calculate the Hall Coefficient
 
            # Hall Coefficient in configuration A
            HC_13_42 = linregress(hall_data[Ti,:,1], hall_data[Ti,:,2])  #Regression on (x,y) = (H_applied, V_measured)
            #Hall Coefficient in configuration B (with source and sense leads swapped)
            HC_24_31 = linregress(hall_data[Ti,:,1], hall_data[Ti,:,4]) #Regression on (x,y) = (H_applied, V_measured)
            
            # Hall coefficient average
            HC_av = linregress(hall_data[Ti,:,1], hall_data[Ti,:,6])
            #print('HC_av',HC_av[0])
            
            ### Step 3: Calculate the charge carrier density and its corresponding error
            
            #calculate the charge carrier density
            cc_density = 1e-6 * np.divide(1, np.multiply(HC_av[0], scipy.constants.e))
            
            ## Error calculation for the Hall coefficient values
            # Taking the error from the linear regression of the average as this is the most accurate value
            Hall_error = HC_av.stderr
            #print('Hall_error',Hall_error)
            # Charge carrier density error is calculated from error propogation
            cc_density_error = 1e-6 *np.divide(Hall_error, np.multiply(HC_av[0]**2, scipy.constants.e))
            #print('cc_density:',cc_density,'cc_density_error',cc_density_error)
            
            #### Step 4: Calculate the mobility and its corresponding error
            
            #extract the resistivity and its error
            resitivity = ppms.res_data[Ti,int(ctf[5]/2),4]
            resistivity_error = ppms.res_data[Ti, int(ctf[5]/2), 5]
        
            # Mobility is calculated from the Hall coefficient and the zero field resistivity (1e4 to convert to cm^2/Vs)
            mobility = 1e4*np.divide(HC_av[0], resitivity)
            
            # Error in the mobility is calculated from error propogation (1e4 to convert to cm^2/Vs)
            # u = Rh/rho -> d(mobility) = sqrt((du/dRh *dRh)^2 + (du/drho *drho)^2)
            mobility_error = 1e4*np.sqrt((Hall_error/resitivity)**2 + ((HC_av[0]*resistivity_error)/resitivity**2)**2)
            
            
            # Average the temperatures over all field values to get a measurement average temperature
            average_temperature = np.sum(tf_av[Ti,:,0], axis=0)/tf_av.shape[1] 
            
            hall_coefficient[Ti,:] = [average_temperature, HC_13_42[0], HC_13_42[2], HC_24_31[0],HC_24_31[2],HC_av[0],HC_av[2], cc_density, cc_density_error, mobility, mobility_error]
            
        # Flatten the hall_data array to a 2D array so it can be put into a df for debugging
        hall_data_flat = np.copy(hall_data).reshape((ctf[4]*ctf[5],7))
        
        # Convert the numpy array to a pandas dataframe for the Hall resistivity
        hall_data_df = pd.DataFrame(hall_data_flat, columns=['Temp (K)', 'Field (T)', 'rho_xy_A(ohm.m)', 'R_squared(I)_A', 'rho_xy_B(ohm.m)','R_squared(I)_B', 'rho_xy_average(ohm.m)'])
        
        # Convert the numpy array to a pandas dataframe for the Hall coefficeint
        hall_coefficient_df = pd.DataFrame(hall_coefficient, columns=['Temp (K)', 'Hallco_A', 'R^2(H)_A', 'Hallco_B','R^2(H)_B', 'Hallco_average','R^2(H)_average', 'Charge Carrier Density (cm^-2)', 'Charge Carrier Density Error (cm^-2)', 'Mobility (cm^2/Vs)', 'Mobility Error (cm^2/Vs)'])
        
        # Store the data in the PPMSData object
        ppms.hall_data = hall_data
        ppms.hall_data_df = hall_data_df
        ppms.hall_coefficient = hall_coefficient
        ppms.hall_coefficient_df = hall_coefficient_df
        
    return PPMS_files
