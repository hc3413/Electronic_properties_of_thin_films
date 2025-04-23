#import all the libraries needed
from import_dep import * # Import Dependencies
from Class_Import import * # Import the data storage class
from Functions_General import * # Import the general functions used to process the data

### ----------- Functions for extraction and analysis of Hall Bar data ----------- ###
 
    

def hallbar_resistivity(
    PPMS_files, #list of PPMSData objects
    print_val = False, # print the initial guess vs calculated resistivity
    filt_kern = 0, # median filter kernel size, if set to 0 then no filter is applied
    filt_sigma = 0, # gaussian filter sigma, if set to 0 then no filter is applied
    threshold = 0 # z-score threshold for outlier detection, if set to 0 then no filter is applied, (lower threshold means more points are considered outliers (2 is typical value))
    ):
    '''Calculates the bulk resistivity (rho_xx) of a thin film sample which has been processed into a Hall bar geometry.
    It first calculates the sheet resistance (R_sheet = R * w/l) and then multiplies by thickness (rho_xx = R_sheet * t).
    The measurements of voltage parallel to the current from the arms on either side are stored in indices 0,1 of the data array.
    Two configurations were used for robustness, to enable a comparison of the resistivity values obtained from the arms on either side of the Hall bar.
    '''
    
    for ppms in PPMS_files:
        
        # Check if the measurement type is 'HallBar' otherwise skip the processing
        if ppms.measurement_type == 'HallBar':
            
            # Extract the required data from the PPMSData object
            film_thickness = ppms.film_thickness
            ctf = ppms.ctf
            tf_av = ppms.tf_av
            data_np = ppms.data_np
            data_np_nd = ppms.data_np_nd
            
            
            # Initialise a list to store the R^2 value for the IV data and check we have a good fit for the linear regression
            R_squared = []
            
            
            # Initialize an empty np aray with indices: (temp_index, field_index, data_colums) 
            # storing each temperature, field, and the corresponding bulk resitivities (rho_xx) for: config A, config B, average of A and B along with the error
            res_data = np.zeros((ctf[4],ctf[5], 6)) # Use 6 columns: Temp, Field, rho_A, rho_B, rho_avg, rho_error

            #Loop over each temperature and field combination, calculating the sheet resistivity using the Van der Pauw method
            for Ti in range(ctf[4]): #for each temperature index
                for Bi in range(ctf[5]): #for each field index
                    
                    ##### Step 1: Using linear regression on the current-voltage data to get resistance R = V/I
                    # R_ij_kl contains: slope(resistance), intercept, R-squared value, p-value, stderr

                    # First pair of Hall bar arms parallel to the current
                    R_A_fit = linregress(data_np_nd[Ti,Bi,:,2,0], data_np_nd[Ti,Bi,:,4,0])    
                    
                    #Second pair of Hall bar arms parallel to the current
                    R_B_fit = linregress(data_np_nd[Ti,Bi,:,2,1], data_np_nd[Ti,Bi,:,4,1])


                    # Append the R-squared value to the list
                    R_squared.extend([R_A_fit.rvalue**2, R_B_fit.rvalue**2]) # Use rvalue**2 for R-squared
                    


                    ##### Step 2: Calculate the sheet resistance (R_sheet) for a hall bar:
                    hb_width = ppms.hb_dimensions[0] * 1e-6 # width of the Hall bar in meters
                    hb_length = ppms.hb_dimensions[1] * 1e-6 # length between the Hall bar arms in meters
                    geometry_factor = hb_width / hb_length # Dimensionless geometry factor (w/l)

                    # Calculate the sheet resistance from R_A_fit
                    # R_sheet = R * (width/length), units: Ohms
                    Rsheet_A = R_A_fit.slope * geometry_factor 

                    # Calculate the sheet resistance from R_B_fit
                    Rsheet_B = R_B_fit.slope * geometry_factor 

                    #### Step 3: Calculate the error in the sheet resistance for each pair of configurations
                    # Error d(R_sheet) = (w/l) * dR, where dR is the standard error from linregress
                    Rsheet_A_error = R_A_fit.stderr * geometry_factor
                    Rsheet_B_error = R_B_fit.stderr * geometry_factor
                    

                    #### Step 4: Calculate the average sheet resistance and the error in the average
                    
                    # Average the two solutions for the final sheet resistance
                    Rsheet_average = (Rsheet_A + Rsheet_B) / 2

                    # Calculate the error in the average sheet resistance by propagating the errors
                    # dRsheet_avg = 0.5 * sqrt(dRsheet_A^2 + dRsheet_B^2)
                    Rsheet_average_error = 0.5*np.sqrt(Rsheet_A_error**2 + Rsheet_B_error**2)
                    
                    #### Step 5: Calculate bulk resistivity (rho) and its error (except if thickness = 1 in which case rho = R_sheet)
                    # rho = R_sheet * thickness (units: Ohm.m)
                    rho_A = Rsheet_A * film_thickness
                    rho_B = Rsheet_B * film_thickness
                    rho_average = Rsheet_average * film_thickness
                    # Final bulk resistivity error d(rho) = d(R_sheet) * thickness (assuming thickness error is negligible)
                    rho_average_error = Rsheet_average_error * film_thickness 
                    
                    #### Step 6: Insert the new row to the np data array using tf_av for temperature and field values that are averaged over all the currents
                    # Store bulk resistivity (rho_xx) values in Ohm.m
                    res_data[Ti,Bi,:] = [tf_av[Ti,Bi,0], tf_av[Ti,Bi,1], rho_A, rho_B, rho_average, rho_average_error]

                # Option to filter the resistivity vs field data
                if filt_kern != 0 or filt_sigma != 0 or threshold != 0:
                    # Filter the average bulk resistivity
                    res_data[Ti,:,4] = filter_data(res_data[Ti,:,4], filt_kern, filt_sigma, threshold)
                
                
            # Flatten the res_data array to a 2D array so it can be put into a df for debugging
            res_data_flat = res_data.reshape((ctf[4]*ctf[5],6))  
            # Convert the numpy array to a pandas dataframe 
            res_data_df = pd.DataFrame(res_data_flat, columns=['Temp (K)', 'Field (T)', 'rho_xx_A (Ohm.m)', 'rho_xx_B (Ohm.m)','rho_xx_average (Ohm.m)', 'rho_error (Ohm.m)'])
            
            # Store the data in the PPMSData object
            ppms.res_data = res_data
            ppms.res_data_df = res_data_df
            ppms.R_squared_res = R_squared
            
            # print(f'res_data({ppms.filename})')
            # print(res_data_df.head())
        
    return PPMS_files


def hallbar_hall(
    PPMS_files, # list of PPMSData objects
    filt_kern = 0, # median filter kernel size, if set to 0 then no filter is applied
    filt_sigma = 0, # gaussian filter sigma, if set to 0 then no filter is applied
    threshold = 0 # z-score threshold for outlier detection, if set to 0 then no filter is applied, (lower threshold means more points are considered outliers (2 is typical value))             
    ):
    '''Calculates the Hall resistivity (rho_xy) at every temperature and field strength.
    Uses linear regression on the current-voltage data to obtain the Hall resistance (R_Hall = V_H / I).
    Calculates Hall resistivity (rho_xy = R_Hall * thickness).
    Calculates Hall coefficient (R_H) from the slope of rho_xy vs B.
    Calculates carrier density (n) and mobility (mu = R_H / rho_xx).
    The Hall resistivity is calculated for each configuration of the source and sense leads and the two are averaged.
    (rho_xy = thickness*V(hall)/I(source))'''
    #ctf = [current_unique, temp_unique, field_unique, current_no, temp_no, field_no]
    
    for ppms in PPMS_files:
        
        # Check if the measurement type is 'HallBar' otherwise skip the processing
        if ppms.measurement_type == 'HallBar':
            
            # Extract the required data from the PPMSData object
            film_thickness = ppms.film_thickness
            ctf = ppms.ctf
            tf_av = ppms.tf_av
            data_np_nd = ppms.data_np_nd
            
            # Initialize an empty np aray to store the temperature, field, hall resistivity (rho_xy): in config A, R^2 for I-V fit A, rho_xy in config B, R^2 for I-V fit B, and Average rho_xy 
            hall_data = np.zeros((ctf[4],ctf[5], 7))
            
            # Initialize an empty np array to store the Temperature, Hall coefficient A, R^2(rho_xy vs B) A, Hall Coefficient B, R^2(rho_xy vs B) B, average Hall coefficient, R^2(rho_xy vs B) av, carrier density, carrier density error, mobility, mobility error
            hall_coefficient = np.zeros((ctf[4], 11))

            #Loop over each temperature using regression on the hall_resistivity-field data to obtain the Hall coefficient at each temperature
            for Ti in range(ctf[4]):
                #Loop over each field using regression on the current-voltage data to obtain the Hall resistance at each field
                for Bi in range(ctf[5]):
                    
                    #### Step 1: Use linear regression on the current-voltage data to obtain the Hall resistance (R_Hall = V_H / I) at each field
                    # R_fit[0] = the slope(Hall resistance), R_fit[1] = intercept, R_fit[2] = R-value (correlation coefficient)     
                    
                    # Hall Resistance in configuration A (using index 2)
                    R_A_fit = linregress(data_np_nd[Ti,Bi,:,2,2], data_np_nd[Ti,Bi,:,4,2]) #Regression on (x,y) = (I_source, V_measured) 
                    # Hall resistance in configuration B (using index 3)
                    R_B_fit = linregress(data_np_nd[Ti,Bi,:,2,3], data_np_nd[Ti,Bi,:,4,3]) #Regression on (x,y) = (I_source, V_measured)
                    
                    # Average the two Hall resistances
                    R_hall_average = (R_A_fit.slope + R_B_fit.slope) / 2
                    
                    # Calculate Hall Resistivity (rho_xy = R_Hall * thickness) unless film thickness = 1 then rho_xy = R_Hall
                    rho_xy_A = R_A_fit.slope * film_thickness
                    rho_xy_B = R_B_fit.slope * film_thickness
                    rho_xy_average = R_hall_average * film_thickness
                    
                    # Insert the new row to the np data array: Temp, Field, rho_xy_A, R^2_A(I-V), rho_xy_B, R^2_B(I-V), rho_xy_average
                    hall_data[Ti,Bi,:] = [tf_av[Ti,Bi,0], tf_av[Ti,Bi,1], rho_xy_A, R_A_fit.rvalue**2, rho_xy_B, R_B_fit.rvalue**2, rho_xy_average]

                # Option to filter Data Hall resistivity vs field
                if filt_kern != 0 or filt_sigma != 0 or threshold != 0:
                    # Filter the average Hall resistivity
                    hall_data[Ti,:,6] = filter_data(hall_data[Ti,:,6], filt_kern, filt_sigma, threshold)
                
                    
                # Step 2: Calculate the Hall Coefficient (R_H) from slope of rho_xy vs B
    
                # Hall Coefficient from rho_xy_A vs B
                HC_A_fit = linregress(hall_data[Ti,:,1], hall_data[Ti,:,2])  #Regression on (x,y) = (B_applied, rho_xy_A)
                # Hall Coefficient from rho_xy_B vs B
                HC_B_fit = linregress(hall_data[Ti,:,1], hall_data[Ti,:,4]) #Regression on (x,y) = (B_applied, rho_xy_B)
                
                # Hall coefficient from average rho_xy vs B
                HC_av_fit = linregress(hall_data[Ti,:,1], hall_data[Ti,:,6]) #Regression on (x,y) = (B_applied, rho_xy_average)
                
                # Extract Hall coefficients (slopes)
                Hall_Coeff_A = HC_A_fit.slope
                Hall_Coeff_B = HC_B_fit.slope
                Hall_Coeff_average = HC_av_fit.slope
                
                ### Step 3: Calculate the charge carrier density (n = 1 / (R_H * e)) and its corresponding error
                
                # Calculate the charge carrier density using the average Hall coefficient
                # Factor 1e-6 converts from m^-3 to cm^-3
                cc_density = 1e-6 / (Hall_Coeff_average * scipy.constants.e) if Hall_Coeff_average != 0 else np.nan
                
                ## Error calculation for the Hall coefficient values
                # Taking the error from the linear regression of the average rho_xy vs B fit
                Hall_Coeff_average_error = HC_av_fit.stderr
                
                # Charge carrier density error is calculated from error propagation: dn = |(dn/dR_H) * dR_H| = |-1/(R_H^2 * e) * dR_H|
                # Factor 1e-6 converts units
                cc_density_error = 1e-6 * np.abs(Hall_Coeff_average_error / ((Hall_Coeff_average**2) * scipy.constants.e)) if Hall_Coeff_average != 0 else np.nan
                
                #### Step 4: Calculate the mobility (mu = R_H / rho_xx) and its corresponding error
                
                # Extract the zero field bulk resistivity (rho_xx) and its error from the corrected res_data
                # Assuming zero field is the middle point of the field sweep
                zero_field_index = int(ctf[5]/2) 
                resistivity_bulk_zero_field = ppms.res_data[Ti, zero_field_index, 4] # This now correctly holds rho_xx_average
                resistivity_bulk_zero_field_error = ppms.res_data[Ti, zero_field_index, 5] # This now correctly holds rho_error
            
                # Mobility is calculated from the average Hall coefficient and the zero field bulk resistivity
                # Factor 1e4 converts from m^2/Vs to cm^2/Vs
                mobility = 1e4 * Hall_Coeff_average / resistivity_bulk_zero_field if resistivity_bulk_zero_field != 0 else np.nan
                
                # Error in the mobility is calculated from error propagation: d(mu) = sqrt( (d(mu)/dR_H * dR_H)^2 + (d(mu)/d(rho) * d(rho))^2 )
                # d(mu)/dR_H = 1/rho_xx
                # d(mu)/d(rho) = -R_H / rho_xx^2
                # Factor 1e4 converts units
                if resistivity_bulk_zero_field != 0:
                    mobility_error = 1e4 * np.sqrt(
                        (Hall_Coeff_average_error / resistivity_bulk_zero_field)**2 + 
                        ((Hall_Coeff_average * resistivity_bulk_zero_field_error) / (resistivity_bulk_zero_field**2))**2
                    )
                else:
                    mobility_error = np.nan
                    
                
                
                # Average the temperatures over all field values to get a measurement average temperature
                average_temperature = np.mean(tf_av[Ti,:,0], axis=0) 
                
                # Store: Temp, R_H_A, R^2(B)_A, R_H_B, R^2(B)_B, R_H_av, R^2(B)_av, n, dn, mu, dmu
                hall_coefficient[Ti,:] = [
                    average_temperature, 
                    Hall_Coeff_A, HC_A_fit.rvalue**2, 
                    Hall_Coeff_B, HC_B_fit.rvalue**2, 
                    Hall_Coeff_average, HC_av_fit.rvalue**2, 
                    cc_density, cc_density_error, 
                    mobility, mobility_error
                ]
                
            # Flatten the hall_data array to a 2D array so it can be put into a df for debugging
            hall_data_flat = np.copy(hall_data).reshape((ctf[4]*ctf[5],7))
            
            # Convert the numpy array to a pandas dataframe for the Hall resistivity
            hall_data_df = pd.DataFrame(hall_data_flat, columns=['Temp (K)', 'Field (T)', 'rho_xy_A (Ohm.m)', 'R^2(I-V)_A', 'rho_xy_B (Ohm.m)','R^2(I-V)_B', 'rho_xy_average (Ohm.m)'])
            
            # Convert the numpy array to a pandas dataframe for the Hall coefficient and derived quantities
            hall_coefficient_df = pd.DataFrame(hall_coefficient, columns=['Temp (K)', 'HallCoeff_A (m3/C)', 'R^2(B)_A', 'HallCoeff_B (m3/C)','R^2(B)_B', 'HallCoeff_Avg (m3/C)','R^2(B)_Avg', 'n (cm^-3)', 'n_Error (cm^-3)', 'Mobility (cm^2/Vs)', 'Mobility_Error (cm^2/Vs)'])
            
            # Store the data in the PPMSData object
            ppms.hall_data = hall_data
            ppms.hall_data_df = hall_data_df
            ppms.hall_coefficient = hall_coefficient
            ppms.hall_coefficient_df = hall_coefficient_df
            
    return PPMS_files




