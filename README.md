# Electronic_properties_of_thin_films
Hall and resitivity measurments to determine the mobility and carrier density in thin films. 
In the main part for Van der Pauw geometry measurments on the PPMS.
Extracts data, combines measurements to remove offsets and plots results giving resistivity and mobility of materials. 

- The program is run in "Analysis" juptyer script with the main functions stored and run from functions.py 
- In general for every np data array created (e.g. "data_np") there is a corresponding dataframe: "data_np_df" so the data can be quickly visualised and checked (if the np data is 3d, then only one of the index is chosen for the df)




Make sure you have activated the virtual environement EPOTFvenv which contains all the required packages to run the code


Contact positions in square VDP config: 2 = top left, 3 = top right, 0 = bottom left, 1 = bottom right

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

    

##### Hall Measurement

- V^{B+}_{02,31} = (V^{B+}_{31}(I^+_{02}) - V^{B+}_{31}(I^-_{02})) / 2
- V^{B-}_{02,31} = (V^{B-}_{31}(I^+_{02}) - V^{B-}_{31}(I^-_{02})) / 2
- V_{AHall} = (V^{B+}_{02,31} - V^{B-}_{02,31}) / 2

- V^{B+}_{13,20} = (V^{B+}_{20}(I^+_{13}) - V^{B+}_{20}(I^-_{13})) / 2
- V^{B-}_{13,20} = (V^{B-}_{20}(I^+_{13}) - V^{B-}_{20}(I^-_{13})) / 2
- V_{BHall} = (V^{B+}_{13,20} - V^{B-}_{13,20}) / 2

- **Net Hall Voltage**
    - V_{Hall} = (V_{AHall} + V_{BHall}) / 2