# Electronic Properties of Thin Films Analysis

This repository contains Python-based analysis tools for electronic transport measurements of thin films, including resistivity, Hall effect, and mobility characterization using both Hall Bar and Van der Pauw geometries.

## Overview

The analysis suite processes data from Physical Property Measurement System (PPMS) measurements to extract fundamental electronic transport properties:
- **Resistivity** (temperature and magnetic field dependent)
- **Hall coefficient** and **charge carrier density**
- **Mobility** (charge carrier mobility)
- **Transport mechanism analysis** (Variable Range Hopping, Small Polaron models)

## Analysis Methods

### Van der Pauw (VDP) Method

The Van der Pauw method enables resistivity and Hall effect measurements on arbitrarily shaped thin film samples using four electrical contacts placed at the sample periphery.

#### VDP Resistivity Analysis (`Functions_VDP.py`)

**Methodology:**
- Uses four measurement configurations (indices 2,3,4,5) grouped into two orthogonal pairs
- Solves the Van der Pauw transcendental equation: exp(-π*R_A/R_s) + exp(-π*R_B/R_s) = 1
- Provides robust measurements through configuration averaging and error propagation

**Key Features:**
- Self-consistent multi-configuration measurements
- Uncertainty quantification through error propagation
- Quality assessment via R² values from linear regression
- Temperature and magnetic field dependent analysis

#### VDP Hall Effect Analysis

**Methodology:**
- Two complementary Hall configurations for independent measurement validation
- Linear regression of Hall resistivity vs magnetic field to extract Hall coefficient
- Single-band model for charge carrier density: n = 1/(R_H × e)
- Mobility calculation: μ = R_H/ρ_xx

**Physical Output:**
- Hall coefficient with sign indicating carrier type (positive=holes, negative=electrons)
- 2D charge carrier density (cm⁻²)
- Mobility (cm²/V·s) with uncertainty quantification

### Hall Bar Method

The Hall Bar geometry uses a dedicated sample design with current contacts at the ends and voltage contacts along the sides for precise longitudinal and transverse measurements.

#### Hall Bar Resistivity Analysis (`Functions_HallBar.py`)

**Dual-Stage Regression Methodology:**
1. **Stage 1 (I-V Analysis):** Linear regression on current-voltage data for each temperature/field point to extract resistance with quality assessment (R²)
2. **Stage 2 (Field Dependence):** Analysis of resistance vs magnetic field to extract zero-field resistivity and magnetoresistance

**Advantages:**
- Eliminates contact resistance effects through four-point measurements
- Separates longitudinal (ρ_xx) and transverse (ρ_xy) components
- Direct geometric factor calculation from sample dimensions

#### Hall Bar Hall Effect Analysis

**Direct Hall Measurement:**
- Uses dedicated Hall voltage contacts perpendicular to current flow
- Linear regression of Hall voltage vs magnetic field
- Integrated analysis with resistivity data for mobility calculation

## Transport Model Fitting (`Analysis_Fitting.ipynb`)

### Variable Range Hopping (VRH)
- Fits ln(ρ) vs T^(-1/4) for 2D Mott VRH
- Extracts localization length and density of states at Fermi level
- Includes R² and RMSE fit quality metrics

### Small Polaron Transport
- Temperature-dependent analysis: ln(ρT^(-S)) vs 1000/T
- Distinguishes adiabatic (S=3/2) vs non-adiabatic (S=1) regimes
- Extracts activation energy and polaron mobility prefactor

### Small Polaron Mobility
- Fits mobility data: μ = μ₀/T^s × exp(-E_σ/(k_b T))
- Determines mobility activation energy and scattering mechanisms
- Comprehensive error propagation and Excel export

## Quality Control and Data Processing

**Measurement Validation:**
- R² values from linear regression assess data linearity
- Multi-configuration comparison for consistency checking
- Statistical outlier detection and filtering options

**Error Analysis:**
- Propagated uncertainties throughout analysis chain
- Monte Carlo error propagation for complex calculations
- Uncertainty quantification in all derived parameters

## Software Architecture

### Core Components
- `Class_Import.py`: Data import and PPMS file handling
- `Functions_VDP.py`: Van der Pauw analysis methods
- `Functions_HallBar.py`: Hall Bar analysis methods  
- `Functions_Fitting.py`: Transport model fitting routines
- `Functions_General.py`: Utility functions and data processing

### Analysis Notebooks
- `Analysis_Basic.ipynb`: Data import and basic plotting
- `Analysis_Fitting.ipynb`: Transport mechanism analysis and fitting
- `Analysis_Publication.ipynb`: Publication-quality figure generation

## Usage

### Environment Setup
```bash
# Activate the virtual environment
source EPOTFvenv/bin/activate

# Launch Jupyter
jupyter notebook
```

### Contact Configuration (VDP Square Geometry)
```
Contact positions: 2 = top left, 3 = top right, 0 = bottom left, 1 = bottom right

Measurement configurations:
- R_{32,10}: Current 3→2, voltage 1→0
- R_{20,31}: Current 2→0, voltage 3→1  
- R_{01,23}: Current 0→1, voltage 2→3
- R_{13,02}: Current 1→3, voltage 0→2
```

#### Detailed Van der Pauw Equations

**Resistivity Measurements:**

- R_{32,10} = (V_{10}(I^+_{32}) - V_{10}(I^-_{32})) / (I^+_{32} - I^-_{32})
    - R_{32,10} = (V_sense(I+)[index 2] - V_sense(I-)[index 2]) / (I_source(I+)[index 2] - I_source(I-)[index 2])

- R_{20,31} = (V_{31}(I^+_{20}) - V_{31}(I^-_{20})) / (I^+_{20} - I^-_{20})
    - R_{20,31} = (V_sense(I+)[index 3] - V_sense(I-)[index 3]) / (I_source(I+)[index 3] - I_source(I-)[index 3])

- rho^A_sheet = (pi * f / ln(2)) * (R_{32,10} + R_{20,31}) / 2
- (q - 1) / (q + 1) = (f * cosh^(-1)(e^(ln(2)/f) / 2)) / ln(2)
- Where f can be taken from tables and q = max(R_{32,10}/R_{20,31}, R_{20,31}/R_{32,10})

- R_{01,23} = (V_{23}(I^+_{01}) - V_{23}(I^-_{01})) / (I^+_{01} - I^-_{01})
    - R_{01,23} = (V_sense(I+)[index 4] - V_sense(I-)[index 4]) / (I_source(I+)[index 4] - I_source(I-)[index 4])

- R_{13,02} = (V_{01}(I^+_{13}) - V_{02}(I^-_{13})) / (I^+_{13} - I^-_{13})
    - R_{20,31} = (V_sense(I+)[index 5] - V_sense(I-)[index 5]) / (I_source(I+)[index 5] - I_source(I-)[index 5])
    
- rho^A_sheet = (pi * f / ln(2)) * (R_{01,23} + R_{13,02}) / 2
- (q - 1) / (q + 1) = (f * cosh^(-1)(e^(ln(2)/f) / 2)) / ln(2)
- Where f can be taken from tables and q = max(R_{01,23}/R_{13,02}, R_{13,02}/R_{01,23})

- **Final resistivity**
    - rho = (rho^A_sheet + rho^B_sheet) / 2

**Hall Measurements:**

- V^{B+}_{02,31} = (V^{B+}_{31}(I^+_{02}) - V^{B+}_{31}(I^-_{02})) / 2
- V^{B-}_{02,31} = (V^{B-}_{31}(I^+_{02}) - V^{B-}_{31}(I^-_{02})) / 2
- V_{AHall} = (V^{B+}_{02,31} - V^{B-}_{02,31}) / 2

- V^{B+}_{13,20} = (V^{B+}_{20}(I^+_{13}) - V^{B+}_{20}(I^-_{13})) / 2
- V^{B-}_{13,20} = (V^{B-}_{20}(I^+_{13}) - V^{B-}_{20}(I^-_{13})) / 2
- V_{BHall} = (V^{B+}_{13,20} - V^{B-}_{13,20}) / 2

- **Net Hall Voltage**
    - V_{Hall} = (V_{AHall} + V_{BHall}) / 2

### Data Analysis Workflow
1. Import PPMS data using `Class_Import.py`
2. Process resistivity using appropriate geometry functions
3. Analyze Hall effect for carrier properties
4. Apply transport model fitting for mechanism identification
5. Generate publication figures with error analysis

## Output Data Format

**Resistivity Data:**
- Temperature (K), Magnetic Field (T)
- ρ_xx (Ω·m) with uncertainty
- Configuration-specific values for validation

**Hall Data:**
- Hall coefficient (m³/C) with uncertainty
- Charge carrier density (cm⁻²)  
- Mobility (cm²/V·s) with uncertainty
- R² values for fit quality assessment

**Transport Fitting:**
- Model parameters with confidence intervals
- Fit quality metrics (R², RMSE)
- Excel export for further analysis

## References

- Van der Pauw, L. J. "A method of measuring specific resistivity and Hall effect of discs of arbitrary shape." Philips Research Reports 13.1 (1958): 1-9.
- Mott, N. F. "Variable range hopping." Philosophical Magazine 19.160 (1969): 835-852.
- Holstein, T. "Studies of polaron motion: Part II. The small polaron." Annals of Physics 8.3 (1959): 325-342.