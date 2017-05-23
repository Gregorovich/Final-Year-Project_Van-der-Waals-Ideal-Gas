# Final-Year-Project_Van-der-Waals-Ideal-Gas

There are 2 MATLAB "programs" here, vdwsolve.m and lv_coexistence.m . lv_coexistence.m calls other scripts to make the code more readable.

vdwsolve.m
    
    MATLAB code for finding the dimensionless van der Waal's equation, internal energy,
    enthalpy and Gibbs free energy
    
    This script rearranges the van der Waal's Equation to remove a and b and make it
    independant of the gas. Also solves for Internal Energy, Enthalpy, and Gibbs free
    energy, and plots the graph of P vs \rho and P vs V
    
    This is done by equating dP/dp = 0 and d2P/dp2 = 0 , and solving.
    
    
lv\_coexistence.m (and Isotherms.m , LV_P_cubic.m , LV_PG_smltns.m )
    
    MATLAB code for creating van der Waal's isotherms at a range of temperatures, creating
    the binodal by solving the van der Waal's cubic, and creating the binodal by solving
    P_g = P_l and G_g = G_l simultaneously. This calls 3 other scripts to plot each graph;
    Isotherms.m , LV_P_cubic.m , and LV_PG_smltns.m.
