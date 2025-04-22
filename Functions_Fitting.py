#import all the libraries needed
from import_dep import *


# Parent fitting function that fits the input data to the function chosen - where various functions are defined below
def fit_function(data, function_type, degree=None):
    
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