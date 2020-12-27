from equ import *
import matplotlib.pyplot as plt

def d_dx(f,x0):
    dx=1e-4
    return (f(x0+dx)-f(x0-dx))/(2*dx)

N=3600# 48*3600
tTotal=100# 48*3600.0
dt=tTotal/N

prim = primario()
estr = (Estruct(prim = prim, dt = dt), False)
secr = (Secr(prim = prim, dt = dt), False)
nucleo = (Nucleo(prim = prim, dt = dt), False)
gv = (GV(prim = prim, dt = dt), False)
loca = (LOCA(prim = prim, dt = dt), False)
sie = (SIE(prim = prim, dt = dt), False)
pot_sys = (pot_cts(prim = prim, dt = dt), True)

prim.sis.append(estr)
prim.sis.append(secr)
prim.sis.append(nucleo)
prim.sis.append(gv)
prim.sis.append(loca)
prim.sis.append(sie)
prim.sis.append(pot_sys)

for i in range(1,N):
    mp_tot = prim.MP_neta
    Q_tot = prim.Q_neta
    P = prim.P[-1]
    dp_i = (Q_tot - prim.A(P)*mp_tot)/(prim.VT*d_dx(prim.B,P)+prim.m[-1]*d_dx(prim.A,P)-prim.VT*1e3)
    prim.P.append(P + dt * dp_i)
    new_P = prim.P[-1]
    prim.T.append(T_sat(new_P))
    prim.m.append(prim.m[-1] + dt * mp_tot)
    new_m = prim.m[-1]
    new_mf = new_m - (prim.VT - new_m * vf(new_P))/vfg(new_P)
    new_mg = (prim.VT - new_m * vf(new_P))/vfg(new_P)
    prim.mf.append(new_mf)
    prim.mg.append(new_mg)
    prim.vol_f.append(new_mf*vf(new_P))
    prim.vol_g.append(new_mg*vf(new_P))
    
    prim.time.append(prim.time[-1] + dt)
    
print(prim.P)