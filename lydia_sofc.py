from scipy.integrate import solve_ivp
import numpy as np
from matplotlib import pyplot as plt
from math import exp, log
import itertools  

# Given Values:
F = 96485e3  # Faraday's constant, C/kmol equivalent charge
R = 8314.5  # Gas constant, J/kmol-K

C_dl_an = 1e10  # F/m3 (in PEMFC example this was F/m2.... need to multiply by thickness?)

"Assumed molar fractions"
X_H_Ni = 0.9
X_Vac_Ni = 0.1
X_OH_YSZ = 0.4
X_O2ion_YSZ = 0.4
X_Vac_YSZ = 0.1
X_H2O_YSZ = 0.1
X_O_LSM = 0.8
X_Vac_LSM = 0.2

"Species standard-state thermo T 800C"
"Anode"
h_H_Ni_o = -9.1e6  # standard-state gibbs energy for H adsorbed on Ni surface (J/kmol)
s_H_Ni_o = 80.2e3  # standard-state gibbs energy for H adsorbed on Ni surface (J/kmol-K)
h_O2ion_YSZ_o = -122.3e6  # standard-state gibbs energy for O2- ion adsorbed on YSZ surface (J/kmol)
s_O2ion_YSZ_o = 128.9e3  # standard-state gibbs energy for O2- ion adsorbed on YSZ surface (J/kmol-K)
h_Vac_Ni_o = 20.0e6  # standard-state gibbs energy for a vacancy on Ni surface (J/kmol)
s_Vac_Ni_o = 33.1e3  # standard-state gibbs energy for a vacancy on Ni surface (J/kmol-K)
h_OH_YSZ_o = 47.8e6  # standard-state gibbs energy for OH adsorbed on YSZ surface (J/kmol)
s_OH_YSZ_o = 165.9e3  # standard-state gibbs energy for OH adsorbed on YSZ surface (J/kmol-K)
h_H2O_YSZ_o = -293.4e6  # standard-state gibbs energy for H2O adsorbed on YSZ surface (J/kmol)
s_H2O_YSZ_o = 177.0e3  # standard-state gibbs energy for H2) adsorbed on YSZ surface (J/kmol-K)



E_act_fwd_an_1 = 90.0e6  # Activation energy (J/kmol)
E_act_fwd_an_2 = 90.0e6  # Activation energy (J/kmol)
A_fwd_an_1 = 1e13  # (kmol-cm-s)
A_fwd_an_2 = 1e12  # (kmol-cm-s)
alpha_fwd_an_1 = 0.5
alpha_fwd_an_2 = 0.5

gamma_surf_an = 1.6e-6
Prod_C_ac_an_1 = (X_Vac_Ni * X_OH_YSZ) / (X_H_Ni * X_O2ion_YSZ)
Prod_C_ac_an_2 = (X_H2O_YSZ) / (X_H_Ni * X_OH_YSZ)

temps = np.array([800+273.15, 900+273.15, 1000+273.15])
i_Far_an = np.zeros_like(temps)

phi_an = 0
phi_elyte_0 = 0.6

i_ext = 1000
SV_0 = np.array([phi_elyte_0 - phi_an])
time_span = np.linspace(0, 10, num=10)
i_ext_array = np.linspace(1e-8,200,num=10)

n_an = -1
dSV_dt = np.zeros_like(SV_0)


def residual(t, SV):
    delta_phi_an = SV[0]

    "temperature dependent"
    g_H_Ni_o = h_H_Ni_o - T * s_H_Ni_o
    g_O2ion_YSZ_o = h_O2ion_YSZ_o - T * s_O2ion_YSZ_o
    g_Vac_Ni_o = h_Vac_Ni_o - T * s_Vac_Ni_o
    g_OH_YSZ_o = h_OH_YSZ_o - T * s_OH_YSZ_o
    g_H2O_YSZ_o = h_H2O_YSZ_o - T * s_H2O_YSZ_o

    Delta_g_circ_an_1 = g_Vac_Ni_o + g_OH_YSZ_o - (g_H_Ni_o + g_O2ion_YSZ_o)
    Delta_g_circ_an_2 = g_H2O_YSZ_o + g_Vac_Ni_o - (g_H_Ni_o + g_OH_YSZ_o)

    k_star_fwd_an_1 = A_fwd_an_1 * exp(-E_act_fwd_an_1 / (R * T))  # DeCaluwe eq 13
    k_star_fwd_an_2 = A_fwd_an_2 * exp(-E_act_fwd_an_2 / (R * T))
    
    Delta_g_an_1 = Delta_g_circ_an_1 + R * T * log(Prod_C_ac_an_1)  # Fuller eq 2.10
    Delta_g_an_2 = Delta_g_circ_an_2 + R * T * log(Prod_C_ac_an_2)
   
    "temperature and delta_phi dependent"
    k_fwd_an_1 = k_star_fwd_an_1 * exp(-alpha_fwd_an_1 * n_an * F * delta_phi_an / (R * T))  # DeCaluwe eq 15
    k_fwd_an_2 = k_star_fwd_an_2 * exp(-alpha_fwd_an_2 * n_an * F * delta_phi_an / (R * T))


    k_rev_an_1 = k_fwd_an_1 * exp(Delta_g_an_1 / (R * T))*exp(n_an*F*delta_phi_an/(R*T))
    k_rev_an_2 = k_fwd_an_2 * exp(Delta_g_an_2 / (R * T))*exp(n_an*F*delta_phi_an/(R*T))

    q_dot_rop_an_1 = (k_fwd_an_1 * (X_H_Ni * X_O2ion_YSZ) - k_rev_an_1 * (X_Vac_Ni * X_OH_YSZ))*gamma_surf_an**2  # DeCaluwe eq 17
    q_dot_rop_an_2 = (k_fwd_an_2 * (X_H_Ni * X_OH_YSZ) - k_rev_an_2 * X_H2O_YSZ*X_Vac_YSZ)*gamma_surf_an**2

    i_Far_an = F * (n_an * q_dot_rop_an_1 + n_an * q_dot_rop_an_2)
    i_dl_an = i_ext - i_Far_an
    dSV_dt[0] = -i_dl_an / C_dl_an

    return dSV_dt

### Case 1: vary temp

for T in temps:
    solution = solve_ivp(residual, time_span, SV_0, rtol=1e-4, atol=1e-6, method='BDF')
    SV_0 = solution.y[:,-1]
    f1 = solution.t
    f2 = np.transpose(solution.y[0])
    plt.plot(f1, f2)
plt.xlabel('Time (s)', fontsize=14)
plt.ylabel('Electric Potential (V)', fontsize=14)
plt.legend(['Anode, T=800 C', 'Anode, T=900 C', 'Anode, T=1000 C',])
plt.title('Cathode Double Layer Electric Potential at Various SOFC Operating Temperatures')
plt.show()

## Case 2: Vary external current
f2_app=[]
for i_v in i_ext_array:
    T = 1073
    i_ext=i_v
    solution = solve_ivp(residual, time_span, SV_0, rtol=1e-4, atol=1e-6,method ='BDF')
    # SV_0 = solution.y[:,-1]
    f1 = i_ext_array
    f2 = solution.y[0,-1]
    f2_app.append(f2)
    f2_vec=np.stack(f2_app, axis=0)
plt.plot(f1, f2_vec)
plt.xlabel('Current (A/m^2)', fontsize=14)
plt.ylabel('Electric Potential (V)', fontsize=14)
plt.show() # change for current on bottom

## Case 3: Vary assumed molar fractions on the nickel

vac_array=[0.1, 0.2, 0.3, 0.4, 0.5]
h_array=[0.9, 0.8, 0.7, 0.6, 0.5]

for (a, b) in zip(vac_array, h_array):
    X_Vac_Ni=a
    X_H_Ni=b
    g_H_Ni_o = h_H_Ni_o - T * s_H_Ni_o
    g_O2ion_YSZ_o = h_O2ion_YSZ_o - T * s_O2ion_YSZ_o
    g_Vac_Ni_o = h_Vac_Ni_o - T * s_Vac_Ni_o
    g_OH_YSZ_o = h_OH_YSZ_o - T * s_OH_YSZ_o
    g_H2O_YSZ_o = h_H2O_YSZ_o - T * s_H2O_YSZ_o
    Delta_g_circ_an_1 = g_Vac_Ni_o + g_OH_YSZ_o - (g_H_Ni_o + g_O2ion_YSZ_o)
    Delta_g_circ_an_2 = g_H2O_YSZ_o + g_Vac_Ni_o - (g_H_Ni_o + g_OH_YSZ_o)
    Prod_C_ac_an_1 = (X_Vac_Ni * X_OH_YSZ) / (X_H_Ni * X_O2ion_YSZ)
    Prod_C_ac_an_2 = (X_H2O_YSZ) / (X_H_Ni * X_OH_YSZ)
    Delta_g_an_1 = Delta_g_circ_an_1 + R * T * log(Prod_C_ac_an_1)  # Do append to add these to a vector
    Delta_g_an_2 = Delta_g_circ_an_2 + R * T * log(Prod_C_ac_an_2)
    solution = solve_ivp(residual, time_span, SV_0, rtol=1e-4, atol=1e-6,method ='BDF')
    # SV_0 = solution.y[:,-1]
    f1 = solution.t
    f2 = np.transpose(solution.y[0])
    plt.plot(f1, f2)
plt.xlabel('Time (s)', fontsize=14)
plt.ylabel('Electric Potential (V)', fontsize=14)
plt.legend(['0.1/0.9', '0.2/0.8','0.3/0.7','0.4/0.6','0.5/0.5'], loc='upper left')
plt.title('Anode Double Layer Electric Potential at Various Concentrations of Vacancies and Hydrogens on Nickel')
plt.show()


