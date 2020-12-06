
from inputs import *
import numpy as np

# Initialize:
phi_dl_an_0 = phi_an_0 - phi_sep_0
phi_dl_ca_0 = phi_ca_0 - phi_sep_0


capacity_anode = H_an * (capacity_gr*eps_gr*density_gr + capacity_si*eps_si*density_si)
capacity_cathode = capacity_nmc*H_ca*eps_nmc*density_nmc
capacity_area = min(capacity_anode,capacity_cathode)


t_final = charge_frac*3600./C_rate

SV_0 = np.array([phi_dl_an_0, phi_dl_ca_0])

# Load inputs and other parameters into 'pars' class:
class pars:
    T = T

    dPhi_eq_an = dPhi_eq_an
    dPhi_eq_ca = dPhi_eq_ca

    i_o_an = i_o_an
    n_an = n_an
    beta_an = beta_an

    i_o_ca = i_o_ca
    n_ca = n_ca
    beta_ca = beta_ca

    C_dl_an_inv = 1/C_dl_an
    C_dl_ca_inv = 1/C_dl_ca

    i_ext = C_rate*capacity_area

    A_fac_an = (r_gr ** 3 + r_si ** 3) / (r_gr ** 2 + r_si ** 2) / 3 / H_an / (eps_gr + eps_si)  
    A_fac_ca = r_nmc/3/H_ca/eps_nmc


i_ext = C_rate*capacity_area