def CO2E_func(i_ext, SV_0, plot_flag):
    from scipy.integrate import solve_ivp
    import numpy as np
    from matplotlib import pyplot as plt
    from math import exp

    C_dl_an = 1e4   # F/m2
    C_dl_ca = 1e4   # F/m2

    i_o_an = 2.5    # anodic exchange current density [=] A/m2
    i_o_ca = 1      # cathodic exchange current density [=] A/m2
    n_an = 2        # number of electrons transfered for anode reaction
    n_ca = 2        # number of electrons transfered for cathode reaction
    F = 96485       # Faraday's constant [=] C/mol
    beta_ca = 0.5   
    beta_an = 0.5
    R = 8.3145
    T = 298.15

    delta_Phi_eq_an = 0.44
    delta_Phi_eq_ca = 0.11
    
    time_span = np.array([0,100])

    # define a derivative. 
    def residual(t,SV):
        dSV_dt = np.zeros_like(SV)

        # Anode Interface:
        eta_an = SV[0] - delta_Phi_eq_an
        i_Far_an = i_o_an*(exp(-n_an*F*beta_an*eta_an/R/T)
                          - exp(n_an*F*(1-beta_an)*eta_an/R/T))
        i_dl_an = i_ext - i_Far_an
        dSV_dt[0] = -i_dl_an/C_dl_an

        # Cathode Interface:
        eta_ca = SV[1] - delta_Phi_eq_ca
        i_Far_ca = i_o_ca*(exp(-n_ca*F*beta_ca*eta_ca/R/T)
                          - exp(n_ca*F*(1-beta_ca)*eta_ca/R/T))
        i_dl_ca = -i_ext - i_Far_ca
        dSV_dt[1] = -i_dl_ca/C_dl_ca
        return dSV_dt

    solution = solve_ivp(residual,time_span,SV_0,rtol=1e-6, atol=1e-8, method='BDF')

    return solution.y[:,-1]