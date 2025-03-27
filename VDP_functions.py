#import all the libraries needed
from import_dep import *
from data_import import *

def round_to_sf(x, p):
    '''Rounds a number to a specified number of significant figures'''
    x = np.asarray(x)
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    return np.round(x * mags) / mags

def unique_T_values(data_import_np):
    '''extract the number of temperature points used in the measurements
    Finding the unique discrete temperature points accounting for a small variation in the temperature setpoint
    also large temperature ranges - thus can't just use decimal points alone or s.f. to find unique values'''
    # Find the rounded values of temperature used and the number of unique temperature points
    temp_round = np.round(data_import_np[:,0,0],decimals=0)
    temp_unique, unique_indices, temp_unique_counts = np.unique(temp_round, return_index=True, return_counts=True)
    
    # Preserve the order of unique temperatures as they appear in temp_round (i.e. the order they appear in the data thus working for ascending or desceding temp measurements)
    sorted_order = np.argsort(unique_indices)
    temp_unique = temp_unique[sorted_order]
    temp_unique_counts = temp_unique_counts[sorted_order]
    
    # Checking for Erronous temperares that are slightly off from the SET value but are not new temperature setponts
    # Set threshold of 80% of the mean temperature frequency, any temperatures appearing below this frequency are considered an error   
    if temp_unique[temp_unique_counts <= np.mean(temp_unique_counts)*0.80].shape[0] >= 1:
        print('tempshape',temp_unique[temp_unique_counts <= np.mean(temp_unique_counts)*0.8].shape[0])
        print('WARNING Potential Temperature Issues: The erronous temperature points are:',temp_unique[temp_unique_counts <= np.mean(temp_unique_counts)*0.8])
        # Remove the erronous temperature points that don't appear frequently enough thus are probably from small temperature fluctuations
        temp_unique = temp_unique[temp_unique_counts >= np.mean(temp_unique_counts)*0.8]
    
    temp_no = temp_unique.shape[0]
    
    return temp_no, temp_unique

def extract_ctf(PPMS_files,  Reduced_temp = False, Reduced_current = False, ohm_m = False):
    '''extract the number of temperature, field and current points used in the measurements
    Also extract the rounded values of the temperature, field and current used in the measurements
    These rounded values can be displayed to check they are as expected where the true more accurate values are used for the calculations
    tf_av are the true, meausured average temperature and field values used for each set of current measurements
    Reduced_temp = [3,-1] will skip the first 3 temperature points and the last 1 temperature point
    Reduced_current = 2 will skip the first 2 current points and the last 2 current points

    '''
    
    for ppms in PPMS_files:
        data_import_np = ppms.data_import_np
    
        #### Step 1: Extract the rounded Current, Temperature, and Field values used from the data

        # Find the rounded values of current used and the number of unique current points
        current_round = round_to_sf(data_import_np[:,2,0],2) #round to 2 significant figures
        current_unique = np.unique(current_round) #find the unique values
        current_unique[(current_unique > -1e-13) & (current_unique < 1e-13)] = 0 #handles the case where the current is 0 but gives a non-zero value
        current_no = current_unique.shape[0] #find the number of unique values

        # Removing current values that are not wanted to "threshold" the current data
        if Reduced_current != False:
            new_current_no = current_no - 2*Reduced_current
        # Initialize an empty np array to store the sliced data
            data_sliced_np = np.zeros([int(data_import_np.shape[0]*(new_current_no/current_no)),data_import_np.shape[1],data_import_np.shape[2]])
            
            # Clip the starting and ending current values by the index of Reduced_current
            for i in range(Reduced_current, current_no - Reduced_current):
                #Slice the data for each current point and store in correct position in new array (now jumping by new_current_no compared to current_no)
                data_sliced_np[i-Reduced_current::new_current_no,:,:] = data_import_np[i::current_no, :, :] 

            # Update data_import_np with the sliced data
            data_import_np = data_sliced_np
            
            # Update the current variables
            current_no = current_no - 2*Reduced_current
            current_unique = current_unique[Reduced_current:- Reduced_current]
    
        # Extract the unique temperature values and the number of unique temperature points
        [temp_no, temp_unique] = unique_T_values(data_import_np)

        # Find the rounded values of field used and the number of unique field points
        # Can't use the np.unique for field as there are repeats of 0 field so it would miscount the field points used per loop
        field_no = int(data_import_np.shape[0]/temp_no/current_no)
        field_unique = data_import_np[0:current_no*field_no:current_no,1,0]

        if Reduced_temp != False:
            data_import_np = data_import_np[Reduced_temp[0]*current_no*field_no:(temp_no + Reduced_temp[1])*current_no*field_no,:,:]
            
            # Repeat the extraction of the unique temperature values and the number of unique temperature points for the new data
            [temp_no, temp_unique] = unique_T_values(data_import_np)
    

        
        #### Step 2: Re-order and store  the data in multi dimensional np array of data vectors (temp, field, source A, source V, sense V)
        # previously: array has (rows, columns, index) dimensions
        # new: (temperature (ctf[4]), field(ctf[5]), current (ctf[3]), columns (temp, field, source A, source V, sense V), index(6))
        # Parameters for reshaping the data
        columns = 5
        indexes = 6
        # Reshape the data to: (temperature, field, current, columns (temp, field, source A, source V, sense V), index(6))
        data_out_nd = np.reshape(data_import_np, (temp_no, field_no, current_no, columns, indexes))


        #### Step 3: Re-order the extracted field value vector so they are in ascending order from -H max to H max 

        # Case for field values that are originally in the order 0->Bmax,-Bmax->0
        if np.round(field_unique[0],decimals=0) == 0 and np.round(field_unique[-1],decimals=0) == 0:
            print('double: Field values originally in the order 0->Bmax,-Bmax->0')
            # Re-order the extracted field values vector
            field_unique = np.concatenate((field_unique[int(field_no/2):],field_unique[:int(field_no/2)]))
            
            # Re-order the main data array to have the field values in the order -Hmax->Hmax
            data_np_nd_reorder = np.zeros_like(data_out_nd)
            for T in range(temp_no):
                start = 0
                middle = int(field_no/2)
                end = int(field_no)
                data_np_nd_reorder = np.concatenate((data_out_nd[:,middle:end,:,:,:],data_out_nd[:,start:middle,:,:,:]), axis = 1)
            
            # Store the reordered data in the data_out variable
            data_out = data_np_nd_reorder
            
            
        # Case for field values that are originally in the order 0,-Hmax->Hmax
        elif np.round(field_unique[0],decimals=0) == 0 and np.round(field_unique[int(field_no/2)], decimals=0) == 0:
            print('single: Field values originally in the order 0,-Bmax->0->Bmax')

            # Re-order the extracted field values vector
            field_unique = np.concatenate((field_unique[1:int(field_no/2)],np.array([field_unique[0]]),field_unique[int(field_no/2):]))
            
            # Re-order the main data array to have the field values in the order -Hmax->Hmax
            data_np_nd_reorder = np.zeros_like(data_out_nd)
            for T in range(temp_no):
                start = 0
                middle = int(field_no/2)
                end = int(field_no)
                data_np_nd_reorder = np.concatenate((data_out_nd[:,1:middle+1,:,:,:], data_out_nd[:,0:1,:,:,:], data_out_nd[:,middle+1:end,:,:,:]), axis = 1)
            
            # Store the reordered data in the data_out variable
            data_out = data_np_nd_reorder
            
        else:
            print('Warning: no recognised field order therefore no reordering of the field values was applied')
            print('Field values:',field_unique[0],field_unique[-1],field_unique[int(field_no/2)])
            data_out = data_out_nd
      
            
        #### Step 4: Store the current, temperature and field data in an array to output from the function for general use
        ctf = [current_unique, temp_unique, field_unique, current_no, temp_no, field_no]
        print('For file:',ppms.filename)
        print(ctf[3],'Currents (uA):',ctf[0]/1e-6)  
        print(ctf[4],'Temperatures (K):',ctf[1])
        print(ctf[5],'Fields (kOe):',np.round(ctf[2]*10,decimals=0))
        print('Is this correct?')
        
            
        #### Step 5: Extract the temperature and field values averaged over all the index values and currents to give a better measurement average
        
        # Temperatures and Field values only - dimensions: (temp_no, field_no, current_no, columns, indexes)
        temp_field_vals = np.copy(data_out[:,:,:,:2,:])
        # Sum over currents - startingdimensions: (temp_no, field_no, current_no, columns, indexes)
        current_averaged = np.sum(temp_field_vals, axis=2)/np.shape(temp_field_vals)[2]
        # Sum over indexes - starting dimensions: (temp_no, field_no, columns, indexes)
        tf_av = np.sum(current_averaged, axis=3)/np.shape(current_averaged)[3]
        # Final dimensions for tf_av: (temp_no, field_no, columns (temp, field))


        #### Step 6: Convert the numpy array to a pandas dataframe for checking values are correct and return the data
        
        # Flatten the data_out array to a 2D array so it can be put into a df for debugging
        data_out_flat = np.copy(data_out).reshape((ctf[4]*ctf[5]*ctf[3],5,6))
        data_out_df = pd.DataFrame(data_out_flat[:,:,2], columns=['Temp (K)', 'Field (T)', 'Source (A)', 'Source (V)', 'Sense (V)'])    
        
        
        
        #### Step 7: Store the data in the PPMSData object

        ppms.data_np = data_out_flat
        ppms.data_np_nd = data_out
        ppms.data_df = data_out_df
        ppms.ctf = ctf
        ppms.tf_av = tf_av
        
        #### Step 8: Scaling factor to output to plots for sheet resistance or resistivity
        # By default all data is handled in ohm-m
        if PPMS_files[0].film_thickness == 1: # In case of 2DEG, you need sheet resistance 
            # No scaling for sheet resistance
            unit_scale = 1
            
        elif ohm_m == 1: #specific case to force resistivity to be in ohm-m
            unit_scale = 1
            
        else:
            # Convert resistivity from ohm-m to micro-ohm cm for better readability
            unit_scale = 1e8
    
    return PPMS_files, unit_scale

    

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



def magnetoresistance(PPMS_files, exclude_res = False):
    '''Calculates the magnetoresistance of a thin film sample using the Van der Pauw method.
    Outputs the magnetoresistance values for each temperature and field strength.
    Third index stores the magnetoresistance value for rho_A, rho_B, and the average of the two: rho_film to look at any ayssmetries with direction'''
    
    # Calculate the resistivity of the thin film sample using the Van der Pauw method
    # Option to exclude if the data is being filtered so you don't overwrite the filtered data
    if exclude_res == False:
        PPMS_files = vdp_resistivity(PPMS_files)
        
    for ppms in PPMS_files:
        # Extract the required data from the PPMSData object
        ctf = ppms.ctf
        res_data = ppms.res_data 
        
        # Initialize an empty array to store the magnetoresistance data: (temperature, field, data_colums)
        # data_colums: (rho_A, rho_B, and the average of the two: rho_film) to look at any ayssmetries with direction
        mag_res = np.zeros((ctf[4], ctf[5],3))
        
        # Loop over the temperature and field data extracting the magnetoresistance values
        for Ti in range(ctf[4]):
            for Bi in range(ctf[5]):
                # Magnetoresistance = 100*(R(H) - R(0)) / R(0)
                mag_res[Ti,Bi,:] = 100*(res_data[Ti,Bi,2:5]-res_data[Ti,int((ctf[5]/2)-1),2:5])/res_data[Ti,int((ctf[5]/2)-1),2:5]
                # Using the first zero for both field orders which should have come from the later, more stable measurement
    
        # Store the magnetoresistance data in the PPMSData object
        ppms.mag_res = mag_res
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
            resitivity = ppms.res_data[Ti,int((ctf[5]/2)-1),4]
            resistivity_error = ppms.res_data[Ti, int((ctf[5]/2)-1), 5]
        
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

def update_plot_string(PPMS_files):
    '''Update the plot_str variable for each file in the tuple.'''
    for ppms in PPMS_files:
        new_plot_str = input(f"Enter new plot string for {ppms.filename} (current: {ppms.plot_str}): ")
        ppms.plot_str = new_plot_str if new_plot_str else ppms.plot_str #ppms.material + ' - ' + 
        print(f"New plot string for {ppms.filename}: {ppms.plot_str}")
    return PPMS_files


def filter_data(raw_data, filt_kern = 0, filt_sigma = 0, threshold = 0):
    '''Filter the data using:
    Median filter 
    Gaussian smoothing filter 
    Z-score outlier detection
    
    filter is applied within the loops of the data processing functions so only applies 
    to a set of IV meausurements at different fields but a single temperature
    
    '''
    
    # Median Filter 
    if filt_kern != 0:
        dat_filtered = scipy.signal.medfilt(raw_data, kernel_size=filt_kern)
        return dat_filtered
        
    # Gaussian Smoothing Filter
    if filt_sigma != 0:
        dat_smoothed = scipy.ndimage.gaussian_filter1d(raw_data, sigma=filt_sigma)
        return dat_smoothed
    
    # Z-score outlier detection 
    if threshold != 0:
        # Calculate the Z-scores of the data
        z_scores = zscore(raw_data)
        
        # Identify outliers
        outliers = np.abs(z_scores) > threshold
        
        # Replace outliers with interpolated values
        dat_filtered = raw_data.copy()
        dat_filtered[outliers] = np.interp(np.flatnonzero(outliers), np.flatnonzero(~outliers), raw_data[~outliers])    
        return dat_filtered

    