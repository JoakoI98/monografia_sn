from equ import *
from matplotlib.pyplot import figure as Figure


def d_dx(f,x0):
    dx=1e-4
    return (f(x0+dx)-f(x0-dx))/(2*dx)

N=36000# 48*3600
tTotal=36000# 48*3600.0
dt=tTotal/N

prim = primario()
estr = (Estruct(prim = prim, dt = dt), True)
secr = (Secr(prim = prim, dt = dt), False)
nucleo = (Nucleo(prim = prim, dt = dt, Q_th_nom=1.05*100E3), True)
gv = (GV(prim = prim, dt = dt), True)
loca = (LOCA(prim = prim, dt = dt), False)
sie = (SIE(prim = prim, dt = dt), False)
pot_sys = (pot_cts(pot = 0, prim = prim, dt = dt), False)

prim.sis.append(estr)
prim.sis.append(secr)
prim.sis.append(nucleo)
prim.sis.append(gv)
prim.sis.append(loca)
prim.sis.append(sie)
prim.sis.append(pot_sys)


prim.debug = False
for el, f in prim.sis:
    el.debug = False
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

t = prim.time
estr, f= estr
nucleo, f= nucleo
gv, f= gv



fig = Figure(figsize=(12,20))
axes = [fig.add_subplot(2,2,i) for i in range(1,5)]
axes[0].plot(t, prim.P)
axes[1].plot(t, prim.T)
axes[2].plot(t, estr.T)
axes[3].plot(t, nucleo.T)

axes[0].set_title("Presion en el primario")
axes[1].set_title("Temperatura en el primario")
axes[2].set_title("Temperatura en estructuras")
axes[3].set_title("Temperatura en el combustible")



for e in axes:
    e.grid(axis = 'both')
    e.set_ylabel("Temperatura [°C]")
    e.set_xlabel("Tiempo [s]")

axes[0].set_ylabel("Presion [MPa]")

fig.show()

print(f"P prim = {prim.P[-1]}   Tprim = {prim.T[-1]}   Testr = {estr.T[-1]}    Tcomb = {nucleo.T[-1]}")
input()
#print(prim.Q_sis)