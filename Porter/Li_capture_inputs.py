# Inputs:

C_rate = 0.1 # How many charges per hour?

T = 298 #K

#Silver anode properties (!! look up capacity of silver, C_dl, verify r_p and H_an are ok)
r_p_an = 4e-6 #m
phi_an_0 = 0 #V, reference condition
C_dl_an = 1e4 #F/m2
i_o_an = 4.0  #A/m2
n_an = -1
beta_an = 0.5
H_an = 30e-6  #m
density_silver = 10490 #kg/m3
capacity_silver = 350 #Ah/kg
eps_silver = .65
dPhi_eq_an = -1.6 #standard potential for Ag, from Appendix A

phi_sep_0 = 1.8  #V

r_p_ca = 0.3e-6 #m
phi_ca_0 = 4.6  #V
C_dl_ca = 1e4 #F/m2
i_o_ca = 100 #A/m2
n_ca = -1
beta_ca = 0.5
H_ca = 50e-6  #m
density_FP = 3056  #kg/m3
capacity_FP = 175  #Ah/kg
eps_FP = 0.65
dPhi_eq_ca = 2.6

# How deep do we want to charge/discharge?
charge_frac = 0.9