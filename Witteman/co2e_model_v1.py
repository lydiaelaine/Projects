from co2e_function_v1_0 import CO2E_func
from scipy.integrate import solve_ivp
import numpy as np
from matplotlib import pyplot as plt
from math import exp

i_array = np.linspace(0, 20000, 100)
V_cell = np.zeros_like(i_array)

phi_an_0 = 0
phi_elyte_an_0 = 1.5
phi_elyte_ca_0 = 0.6
phi_ca = 0

delta_Phi_eq_an = 0.44
delta_Phi_eq_ca = 0.11
delta_Phi_j = 0.828

SV_0 = np.array([phi_an_0 - phi_elyte_an_0 + delta_Phi_j, phi_elyte_ca_0 - phi_ca])
#   SV_0 = np.array([2.0, 0.8])

SV = np.zeros((SV_0.size,i_array.size))
for j,i_ext in enumerate(i_array):
    plot = 0
    SV[:,j] = CO2E_func(i_ext, SV_0, plot)
    SV_0 = SV[:,j]

V_cell = SV[0,:] + SV[1,:]

fig, ax1 = plt.subplots()
ax1.plot(i_array/1e4,V_cell) 
plt.xlabel('Current Density (A/cm$^2$)',fontsize=14)
plt.ylabel('Applied Potential (V)',fontsize=14)
plt.show()