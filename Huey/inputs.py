# Inputs
C_rate = 0.1 # How many charges per hour?

T = 298 #K

r_gr = 30E-6 #m, from CAMP
r_si = 150E-9 #m, from CAMP 
phi_an_0 = 0 #V
C_dl_an = 1e4 #F/m2
i_o_an = 4.0  #A/m2
n_an = -1
beta_an = 0.5
H_an = 30e-6  #m
density_gr = 2260 #kg/m3
capacity_gr = 350 #Ah/kg
eps_gr = .70 #can go back and actually caluculate this from CAMP data, but i didn't yet
density_si = 2330 #kg/m3
capacity_si = 4200 #Ah/kg
eps_si = .20
dPhi_eq_an = -1.6

phi_sep_0 = 1.8  #V

r_nmc = 0.3e-6 #m
phi_ca_0 = 4.6  #V
C_dl_ca = 1e4 #F/m2
i_o_ca = 100 #A/m2
n_ca = -1
beta_ca = 0.5
H_ca = 50e-6  #m
density_nmc = 2200  #kg/m3, from Targray NMC 811
capacity_nmc = 185  #Ah/kg
eps_nmc = 0.65
dPhi_eq_ca = 2.6

# How deep do we want to charge/discharge?
charge_frac = 0.9