import numpy as np

class prim:
    p_n = 12.25
    pot = 100
    vol_t = 55.9
    vol_a = 13.31
    vol_v = 10.4
    m_rpv = 143
    m_irpv = 53
    
class sie:
    vol_a = 50
    vol_n2 = 19
    T_i = 45
    k_d = 2.47E4
    p = 2
    p_r = 0.5

class comb:
    r1 = 3.8
    r2 = 3.875
    r3 = 4.5
    N_b = 6583
    L_a = 1.4
    r_cp_uo2 = 3E6
    r_cp_v = 2E6

def q_dec(t, q0):
    if t < 9:
        val = 0.15*q0*np.exp(-0.1*t)
    elif t < 180:
        val = 0.09192*q0*t**-0.181
    else:
        val = 0.156*q0*t**-0.283
    return val

class gv:
    p_s = 4.7
    m_l = 600

class water_f:
    
    R = 0.4615
    
    def T_sat(self, p):
        return 178.60744 + 103.36613*np.log10(p) + 29.7073*np.log10(p)**2
    
    def hf(self, p):
        return 759.25648+402.51353*np.log10(p)+144.09837*np.log10(p)**2+114.95627*np.log10(p)**3
    
    def hg(self, p):
        return 2742.65249+16.48065*p-1.66783*p**2
    
    def rhof(self, p):
        return 938.44405-44.48168*p+2.76929*p^2-0.08432*p^3
    
    def rhog(self, p):
        return 2.01394+3.01869*p+0.23865*p^2
    
    def cp(self, T):
        return (32.24+0.1923e-2*(T+273.15)+1.055e-5*(T+273.15)**2-3.595e-9*(T+273.15)**3)/18
    
    def cv(self, T):
        return self.cp(T) - R