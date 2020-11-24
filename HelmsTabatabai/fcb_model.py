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
C_rate = 6 # How many charges per hour? 

T = 298 #K W

charge_frac = 0.9
r_p_an = 4e-6 #m
phi_an_0 = 0 #V 
C_dl_an = 1e4 #F/m2
i_o_an = 4.0  #A/m2
n_an = -1
beta_an = 0.5
H_an = 30e-6  #m
density_graphite = 2260 #kg/m3
capacity_graphite = 350 #Ah/kg
eps_graphite = .65
dPhi_eq_an = -1.6

phi_sep_0 = 1.8  #V W 

r_p_ca = 0.3e-6 #m
phi_ca_0 = 4.6  #V 
C_dl_ca = 1e4 #F/m2
i_o_ca = 100 #A/m2
n_ca = -1
beta_ca = 0.5
H_ca = 50e-6  #m
density_LCO = 2292  #kg/m3
capacity_LCO = 175  #Ah/kg
eps_LCO = 0.65
dPhi_eq_ca = 2.6

F = 96485 #C/mol
R = 8.314 #J/mol*K

#mass transport parameters
mu =  0.0015  #viscosity kg/m/s
n_brugg = -0.5 #Bruggeman coefficient 
d_part = 1e-6 #m, avg diameter of active material particle
C_elyte = 1000 #mol/m3, concetration of electrolye 

"End User Inputs"

# Import necessary modules:
from scipy.integrate import solve_ivp #integration function for ODE system.
import numpy as np
from matplotlib import pyplot as plt
from math import exp

phi_dl_an_0 = phi_an_0 - phi_sep_0
phi_dl_ca_0 = phi_ca_0 - phi_sep_0


capacity_anode = capacity_graphite*H_an*eps_graphite*density_graphite
capacity_cathode = capacity_LCO*H_ca*eps_LCO*density_LCO
capacity_area = min(capacity_anode,capacity_cathode)


t_final = charge_frac*3600./C_rate
i_ext = C_rate*capacity_area

A_fac_an = r_p_an/3/H_an/eps_graphite
A_fac_ca = r_p_ca/3/H_ca/eps_LCO

#Initial solution vector:

SV_0 = np.array( phi_ca_0 - phi_sep_0)
    
def residual(t,SV):
    dSV_dt = np.zeros_like(SV)
    
    eta_an = SV[0] - dPhi_eq_an
    i_Far_an = i_o_an*(exp(-n_an*F*beta_an*eta_an/R/T)
                      - exp(n_an*F*(1-beta_an)*eta_an/R/T))
    i_dl_an = i_ext*A_fac_an - i_Far_an
    dSV_dt[0] = i_dl_an/C_dl_an
    
    
    eta_ca = SV[1] - dPhi_eq_ca
    i_Far_ca = i_o_ca*(exp(-n_ca*F*beta_ca*eta_ca/R/T)
                      - exp(n_ca*F*(1-beta_ca)*eta_ca/R/T))
    i_dl_ca = -i_ext*A_fac_ca - i_Far_ca
    
    
    dSV_dt[1] = i_dl_ca/C_dl_ca
    
    return dSV_dt
    
from scipy.integrate import solve_ivp

SV_0 = np.array([phi_dl_an_0, phi_dl_ca_0])

time_span = np.array([0,t_final])

solution = solve_ivp(residual,time_span,SV_0,rtol=1e-6, atol=1e-8)  

from matplotlib import pyplot as plt
for var in solution.y:
    plt.plot(solution.t,var)
    
plt.legend(['Anode double layer','Cathode double layer'])
plt.ylabel('Potential (V)',fontsize=14)
plt.xlabel('Time(s)',fontsize=14)

# Either directly in this file, or in a separate file that you import, define: 
#   - A residual function called 'residual'
#   - An array 'time_span' which has [0, t_final] where t_final is the total 
#       length of time you want the model to simulate.
#   - An intial solution vector SV_0
