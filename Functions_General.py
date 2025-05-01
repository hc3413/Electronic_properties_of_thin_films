#import all the libraries needed
from import_dep import * # Import Dependencies
from Class_Import import * # Import the data storage class


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

def extract_ctf(
    PPMS_files,  
    Reduced_temp = False,  #Reduced_temp = [3,-1] will skip the first 3 temperature points and the last 1 temperature point
    Reduced_current = False,  #Reduced_current = 2 will skip the first 2 current points and the last 2 current points
    ohm_m = False, # parameter to be passed in to force the unit scaling for plots to be in ohm-m instead of u-ohm cm
    single_zero_field = False, # Removes the first zero field point from the data thus giving a single zero field point for final plotting
    no_of_columns: int = 5  # By default (excluding the index colum which we have removed) there are 5 columns: temp, field, source A, source V, sense V
    ):
    '''extract the number of temperature, field and current points used in the measurements
    Also extract the rounded values of the temperature, field and current used in the measurements
    These rounded values can be displayed to check they are as expected where the true more accurate values are used for the calculations
    tf_av are the true, measured average temperature and field values used for each set of current measurements
    Reshapes the data into a multi-dimensional numpy array: (temperature, field, current, columns, index)
    '''
    
    for ppms in PPMS_files:
        data_import_np = ppms.data_import_np
        measurement_type = ppms.measurement_type # Get measurement type from the object
    
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
            # Ensure slicing indices are valid
            start_idx = Reduced_temp[0] * current_no * field_no
            end_idx = (temp_no + Reduced_temp[1]) * current_no * field_no
            
            # Adjust end_idx if it's negative or exceeds bounds
            if Reduced_temp[1] < 0:
                 end_idx = data_import_np.shape[0] + Reduced_temp[1] * current_no * field_no
            
            # Ensure start_idx is not greater than end_idx
            if start_idx >= end_idx or start_idx < 0 or end_idx > data_import_np.shape[0]:
                 print(f"Warning: Invalid Reduced_temp indices {Reduced_temp} for file {ppms.filename}. Skipping temperature reduction.")
            else:
                 data_import_np = data_import_np[start_idx:end_idx,:,:]
                 # Repeat the extraction of the unique temperature values and the number of unique temperature points for the new data
                 [temp_no, temp_unique] = unique_T_values(data_import_np)
    

        
        #### Step 2: Re-order and store the data in multi dimensional np array of data vectors (temp, field, source A, source V, sense V)
        # previously: array has (rows, columns, index) dimensions
        # new: (temperature (ctf[4]), field(ctf[5]), current (ctf[3]), columns (temp, field, source A, source V, sense V), index(number of measurement configurations))
        # Parameters for reshaping the data
        columns = no_of_columns # Should be 5: Temp, Field, Source I, Source V, Sense V
        
        # Determine the number of indices based on measurement type stored in the object
        if measurement_type == 'VDP':
            indexes = 6 
        elif measurement_type == 'HallBar':
            indexes = 4
        else:
             # Should not happen if import is correct, but good practice to check
             raise ValueError(f"Unknown measurement type '{measurement_type}' in PPMS object for {ppms.filename}")

        # Check if the total number of data points matches the expected shape
        expected_points = temp_no * field_no * current_no
        actual_points = data_import_np.shape[0]
        
        if expected_points != actual_points:
             print(f"Warning: Mismatch in expected ({expected_points}) vs actual ({actual_points}) data points for {ppms.filename}. Check data integrity or ctf calculation.")
             # Optionally, skip this file or try to proceed cautiously
             # continue 
        
        # Reshape the data to: (temperature, field, current, columns, index)
        try:
            data_out_nd = np.reshape(data_import_np, (temp_no, field_no, current_no, columns, indexes))
        except ValueError as e:
             print(f"Error reshaping data for {ppms.filename}: {e}")
             print(f"Expected shape: ({temp_no}, {field_no}, {current_no}, {columns}, {indexes})")
             print(f"Input array size: {data_import_np.size}")
             # Handle error, e.g., skip this file
             continue


        #### Step 3: Re-order the extracted field value vector so they are in ascending order from -H max to H max 
        # ... (Field reordering logic remains the same) ...
        # Case for field values that are originally in the order 0->Bmax,-Bmax->0
        
        # Don't reorder the field values if the rotator is used
        if ppms.rotator == True:
            data_out = data_out_nd # Use original data if no reordering applied
        
        # Case for 'double' field falues with order 0 ->Bmax,-Bmax->0 
        elif np.round(field_unique[0],decimals=0) == 0 and np.round(field_unique[-1],decimals=0) == 0 and field_no > 1:
            print(f'{ppms.filename}: Field values originally in the order 0->Bmax,-Bmax->0')
            
            middle = int(field_no/2)
            # Re-order the extracted field values vector
            if single_zero_field == True:
                # Keep only the final zero field point
                field_unique_reordered = np.concatenate((field_unique[middle:],field_unique[1:middle]))
                field_no_reordered = field_no - 1 # Update field number
            else:
                # Keep both zero field points
                field_unique_reordered = np.concatenate((field_unique[middle:],field_unique[:middle]))
                field_no_reordered = field_no # Field number stays the same
            
            # Re-order the main data array to have the field values in the order -Hmax->Hmax
            data_out_list = []
            for T in range(temp_no):
                start = 0
                end = int(field_no)
                
                if single_zero_field == True:
                    # Concatenate (-Bmax -> 0_end) and (0_start+1 -> Bmax)
                    reordered_slice = np.concatenate((data_out_nd[T,middle:end,:,:,:],data_out_nd[T,start+1:middle,:,:,:]), axis = 0)
                else:
                     # Concatenate (-Bmax -> 0_end) and (0_start -> Bmax)
                    reordered_slice = np.concatenate((data_out_nd[T,middle:end,:,:,:],data_out_nd[T,start:middle,:,:,:]), axis = 0)
                data_out_list.append(reordered_slice)
            
            # Store the reordered data back into a single numpy array
            data_out = np.stack(data_out_list, axis=0)
            field_unique = field_unique_reordered
            field_no = field_no_reordered # Update field_no for ctf
            
            
        # Case for field values that are originally in the order 0,-Hmax->Hmax->0
        elif np.round(field_unique[0],decimals=0) == 0 and np.round(field_unique[int(field_no/2)], decimals=0) == 0 and field_no > 2:
            print(f'{ppms.filename}: Field values originally in the order 0,-Bmax->0->Bmax')

            middle = int(field_no/2)
            # Re-order the extracted field values vector
            if single_zero_field == True:
                # Keep only the middle zero field point: (-Bmax -> 0_mid -> Bmax)
                field_unique_reordered = np.concatenate((field_unique[1:middle+1],field_unique[middle+1:]))
                field_no_reordered = field_no - 1 # Update field number
            else:
                # Keep both zero field points: (-Bmax -> 0_start -> 0_mid -> Bmax)
                field_unique_reordered = np.concatenate((field_unique[1:middle],np.array([field_unique[0]]),field_unique[middle:]))
                field_no_reordered = field_no # Field number stays the same
            
            # Re-order the main data array to have the field values in the order -Hmax->Hmax
            data_out_list = []
            for T in range(temp_no):
                start = 0
                end = int(field_no)
                
                if single_zero_field == True:
                    # Concatenate (-Bmax -> 0_mid) and (0_mid+1 -> Bmax)
                    reordered_slice = np.concatenate((data_out_nd[T,1:middle+1,:,:,:], data_out_nd[T,middle+1:end,:,:,:]), axis = 0)
                else:
                    # Concatenate (-Bmax -> 0_start) and (0_mid -> Bmax)
                    reordered_slice = np.concatenate((data_out_nd[T,1:middle,:,:,:], data_out_nd[T,0:1,:,:,:], data_out_nd[T,middle:end,:,:,:]), axis = 0)
                data_out_list.append(reordered_slice)

            # Store the reordered data back into a single numpy array
            data_out = np.stack(data_out_list, axis=0)
            field_unique = field_unique_reordered
            field_no = field_no_reordered # Update field_no for ctf
            
        else:
            print(f'Warning for {ppms.filename}: no recognised field order for reordering was applied.')
            print(f'Field values start/mid/end: {field_unique[0] if field_no > 0 else "N/A"}, {field_unique[int(field_no/2)] if field_no > 1 else "N/A"}, {field_unique[-1] if field_no > 0 else "N/A"}')
            data_out = data_out_nd # Use original data if no reordering applied
      
            
        #### Step 4: Store the current, temperature and field data in an array to output from the function for general use
        ctf = [current_unique, temp_unique, field_unique, current_no, temp_no, field_no] # Use updated field_no
        print('For file:',ppms.filename)
        print(ctf[3],'Currents (uA):',ctf[0]/1e-6)  
        print(ctf[4],'Temperatures (K):',ctf[1])
        print(ctf[5],'Fields (kOe):',np.round(ctf[2]*10,decimals=0))
        print('Is this correct?')
        
            
        #### Step 5: Extract the temperature and field values averaged over all the index values and currents to give a better measurement average
        
        # Temperatures and Field values only - dimensions: (temp_no, field_no, current_no, columns, indexes)
        temp_field_vals = np.copy(data_out[:,:,:,:2,:]) # Use data_out which might be reordered
        # Average over currents - starting dimensions: (temp_no, field_no, current_no, columns, indexes)
        current_averaged = np.mean(temp_field_vals, axis=2)
        # Average over indexes - starting dimensions: (temp_no, field_no, columns, indexes)
        tf_av = np.mean(current_averaged, axis=3)
        # Final dimensions for tf_av: (temp_no, field_no, columns (temp, field))


        #### Step 6: Convert the numpy array to a pandas dataframe for checking values are correct and return the data
        
        # Flatten the data_out array to a 2D array so it can be put into a df for debugging
        # Use updated field_no and indexes
        data_out_flat = np.copy(data_out).reshape((ctf[4]*ctf[5]*ctf[3], columns, indexes)) 
        # Example: Select one index (e.g., index 0) for the DataFrame view
        df_index_to_view = 0 
        if df_index_to_view >= indexes:
             df_index_to_view = 0 # Default to index 0 if chosen index is out of bounds
        data_out_df = pd.DataFrame(data_out_flat[:,:,df_index_to_view], columns=['Temp (K)', 'Field (T)', 'Source (A)', 'Source (V)', 'Sense (V)'])    
        
        
        
        #### Step 7: Store the data in the PPMSData object

        ppms.data_np = data_out_flat # Store the flattened 3D array (points, columns, index)
        ppms.data_np_nd = data_out # Store the 5D array (temp, field, current, columns, index)
        ppms.data_df = data_out_df # Store the example DataFrame view
        ppms.ctf = ctf # Store the characteristics [currents, temps, fields, counts]
        ppms.tf_av = tf_av # Store the averaged temp/field array
        
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

    


def magnetoresistance(PPMS_files, exclude_res=False):
    '''
    Calculates the magnetoresistance of a thin film sample using the Van der Pauw method.
    Outputs the magnetoresistance values for each temperature and field strength.
    Now includes an extra column for the magnetoresistance error of the average (ρ_film).
    '''
        
    for ppms in PPMS_files:
        ctf = ppms.ctf
        res_data = ppms.res_data 
       
        # Shape (ctf[4], ctf[5], 4) (temperature, field, columns)
        mag_res = np.zeros((ctf[4], ctf[5], 5)) # Columns: MR_A, MR_B, MR_avg, MR_avg_error, MR_fitted
        
        zero_field_index = int(ctf[5] / 2)  # Middle (zero) index

        for Ti in range(ctf[4]):
            for Bi in range(ctf[5]):
                # R0 for average & its error (at zero field)
                r0 = res_data[Ti, zero_field_index, 4]
                e0 = res_data[Ti, zero_field_index, 5]
                
                # R(H) for average & its error
                rh = res_data[Ti, Bi, 4]
                eh = res_data[Ti, Bi, 5]
                
                # MR for ρ_A, ρ_B, and ρ_average
                # A
                mag_res[Ti, Bi, 0] = 100 * (
                    (res_data[Ti, Bi, 2] - res_data[Ti, zero_field_index, 2])
                    / res_data[Ti, zero_field_index, 2]
                )
                # B
                mag_res[Ti, Bi, 1] = 100 * (
                    (res_data[Ti, Bi, 3] - res_data[Ti, zero_field_index, 3])
                    / res_data[Ti, zero_field_index, 3]
                )
                # Average
                mag_res[Ti, Bi, 2] = 100 * (rh - r0) / r0

                # Error propagation for MR_avg:
                # MR(H) = 100 * (rh - r0) / r0
                # partial wrt rh = 100 / r0
                # partial wrt r0 = -100 * rh / r0^2
                dMR_drh = 100.0 / r0
                dMR_dr0 = -100.0 * rh / (r0**2)
                mr_error = np.sqrt((dMR_drh * eh)**2 + (dMR_dr0 * e0)**2)

                mag_res[Ti, Bi, 3] = mr_error

        ppms.mag_res = mag_res

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

