#import all the libraries needed
from import_dep import *


# Parent fitting function that fits the input data to the function chosen - where various functions are defined below
def fit_function_temp(
    PPMS_files, #list of PPMS data objects to be fitted
    function_type: str, #model function to fit the data to
    degree: int = None, #if model = 'polynomial', degree of the polynomial to fit
    T_range: tuple = None, #tuple of (T_min, T_max) to limit the fitting range
    B_range: tuple = None, #tuple of (B_min, B_max) to limit the fitting range
    ):
    '''
    Fit the input data vs temperature to the specified function type.
    
    This function will call the appropriate function for the model based on the function_type.
    
    Resistivity Models:
    1. 'polynomial'
    2. 'bandgap'
    3. 'polaron'
    4. 'variable_range_hopping'
    
    Charge Carrier Density Models:
    1. 'polynomial'
    2. 'bandgap'
    3. 'polaron'
    
    Mobility Models:
    1. 'polynomial'
    2. 'impurity' -    
    
    '''
    
    for ppms in PPMS_files:
        
        ### Step 1: Get measured data from the PPMS object
        
        res_data = ppms.res_data # (T_index, Bindex,  [Temp, Field, rho_A, rho_B, rho_avg, rho_error, rho_fit])
        hall_coefficient = ppms.hall_coefficient #(T_index, [T, 'Rh_A', 'R^2_A', 'Rh_B','R^2_B', 'Rh_av','R^2_av', 'n', 'n_error', 'u', 'u_error', 'n_fitted', 'u_fitted'])
        
        #mag_res = ppms.mag_res # Columns: MR_A, MR_B, MR_avg, MR_avg_error, MR_fitted

        # Step 2: mask the data based on the T_range and B_range
        if T_range is not None:
            T_min, T_max = T_range
            
            # Mask the resistivity data
            res_mask = (res_data[:, :, 0] >= T_min) & (res_data[:, :, 0] <= T_max)
            res_data = res_data[res_mask]
            
            # Mask the Hall data
            hall_mask = (hall_coefficient[:, 0] >= T_min) & (hall_coefficient[:, 0] <= T_max)
            hall_coefficient = hall_coefficient[hall_mask]
            
        if B_range is not None:
            B_min, B_max = B_range
            
            # Mask the resistivity data
            res_mask = (res_data[:, :, 1] >= B_min) & (res_data[:, :, 1] <= B_max)
            res_data = res_data[res_mask]
            
            # Mask the Hall data
            hall_mask = (hall_coefficient[:, 1] >= B_min) & (hall_coefficient[:, 1] <= B_max)
            hall_coefficient = hall_coefficient[hall_mask]
            
            B_fields = res_data[0, :, 1] # extract the B field values from the first row of the resistivity data
        
        
        
        # Step 3: Fit the resistivity data for each B field 
        for Bindx, B in enumerate(B_fields):
            
            temperature = res_data[:, Bindx, 0]
            resistivity = res_data[:, Bindx, 4]
            
            
            if function_type == 'polynomial':
                # Fit resistivity data
                coeffs, fitted_data = fit_polynomial(res_data[:, Bindx, :], degree)
                
                # Format the polynomial coefficients into a readable string
                equation = poly_to_string(coeffs)
                
                # Store the fitted data and equation in the PPMS object
                ppms.res_fit[Bindx] = {
                    'coeffs': coeffs,
                    'fitted_data': fitted_data,
                    'equation': equation
                }
                

    
    
    
    

    
    return
    
    




def fit_polynomial(data, degree):
    """Fit a polynomial of a given degree to the data."""
    x = data[:, 0]
    y = data[:, 1]
    
    # Fit polynomial
    coeffs = np.polyfit(x, y, degree)
    
    # Generate polynomial function
    poly_func = np.poly1d(coeffs)
    
    # Generate x values for plotting
    x_fit = np.linspace(np.min(x), np.max(x), 100)
    y_fit = poly_func(x_fit)
    
    # put fitted data into np.array of same shape as data
    fitted_data = np.column_stack((x_fit, y_fit))
    
    return coeffs, fitted_data
    


# Create a helper function to format polynomial coefficients as a readable equation
def poly_to_string(coeffs):
    """Convert polynomial coefficients to a readable string with proper formatting."""
    powers = range(len(coeffs) - 1, -1, -1)  # Highest power to lowest
    terms = []
    
    for coef, power in zip(coeffs, powers):
        # Skip near-zero coefficients
        if abs(coef) < 1e-10:
            continue
            
        # Format the coefficient with 2 decimal places
        coef_str = f"{abs(coef):.2f}"
        
        # Format the term based on its power
        if power == 0:    # Constant term
            term = coef_str
        elif power == 1:  # Linear term
            term = f"${coef_str}T$"
        else:             # Higher power terms
            term = f"${coef_str}T^{power}$"
            
        # Add + or - sign
        if coef < 0:
            term = f"- {term}"
        elif terms:  # Not the first term and positive
            term = f"+ {term}"
            
        terms.append(term)
    
    equation = " ".join(terms)
    return equation