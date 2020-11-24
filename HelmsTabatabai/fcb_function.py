"""
    This file runs and executes a your model, calculating the cell/device/system properties of interest as a function of time. 

    The code structure at this level is meant to be very simple, and delegates most functions to lower-level modules, which can be updated with new capabilties, over time.  The code:

        1 - Reads inputs and initializes the model
        
        2 - Calls the residual function and integrates over the user-defined    
            time span.
        3 - The simulation then returns the solution vector at the end of the 
            integration, which can be processed as needed to plot or analyze any quantities of interest.
"""

""" 
Model domain:
Cathode- Li metal oxide
Electrolyte- LiPF6 & Ethylene Carbonate:Dimethyl Carbonate 
        Liquid electrolyte with C_LiPF6 in a X_EC:X_DMC mixture 
Anode- Carbon

"""

"User Inputs"

i_ext =   #A
t_final = 500 #s, charging time
T = 298   #K, simulation temperature

"Electric Potentials"
phi_an_0 = 0  #V
phi_elyte_0 = 1.8  #V 
phi_ca_0 = 4.6  #V

C_dl_an = 1e2  #F/m2  capacitance 
C_dl_ca = 1e2  #F/m2

i_o_an = 4.0  #A/m2  exchange current density
i_o_ca = 100  #A/m2  exchange current density  
n_an = 1
n_ca = 1

beta_an = 0.5 #charge transfer symmetry factor
beta_ca = 0.5

F = 96485 #C/mol
R = 8.314 #J/mol*K

#mass transport parameters
mu =  0.0015  #viscosity kg/m/s
n_brugg = -0.5 #Bruggeman coefficient 
d_part = 1e-6 #m, avg diameter of active material particle
C_elyte = 1000 #mol/m3, concetration of electrolye 

"Anode Carbon Properties"
eps_carbon = 0.45  #volume fraction 
rho_carbon = 2260 #kg/m3  density 
capacity_carbon = 350 #Ah/kg

"Cathode LMO Properties"
eps_LMO = 0.45
rho_LMO = 2000 #kg/m3
capacity_LMO = 150 #Ah/kg 


delta_phieq_an = 0.5
delta_phieq_ca = 0.5

"End User Inputs"

# Import necessary modules:
from scipy.integrate import solve_ivp #integration function for ODE system.
import numpy as np

#Initial solution vector"

SV_0 = np.array('[phi_elyte_0 - phi_an_0,' phi_ca_0 - phi_elyte_0)

class pars: 
        time_span = np.array([0,t_final])
        
        T = T
        i_ext = i_ext
        mu = mu
        
        
        #Anode
        delta_phieq_an = delta_phieq_an
        i_o_an = i_o_an
        n_an = n_an
        beta_an = beta_an 
        C_dl_an = C_dl_an 
        eps_carbon = eps_carbon 
        rho_carbon = rho_carbon 
        capacity_carbon = capacity_carbon
        
         #Cathode
        delta_phieq_ca = delta_phieq_ca
        i_o_ca = i_o_ca
        n_ca = n_ca
        beta_ca = beta_ca 
        C_dl_ca = C_dl_ca 
        eps_LMO = eps_LMO 
        rho_LMO = rho_LMO 
        capacity_LMO = capacity_LMO

# Either directly in this file, or in a separate file that you import, define: 
#   - A residual function called 'residual'
#   - An array 'time_span' which has [0, t_final] where t_final is the total 
#       length of time you want the model to simulate.
#   - An intial solution vector SV_0
solution = solve_ivp(residual, time_span, SV_0,rtol=1e-4, atol=1e-6)
