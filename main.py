from equ import *
from matplotlib.pyplot import figure as Figure


def d_dx(f,x0):
    dx=1e-4
    return (f(x0+dx)-f(x0-dx))/(2*dx)

N=3600# 48*3600
tTotal=3600# 48*3600.0
dt=tTotal/N

prim = primario()
estr = (Estruct(prim = prim, dt = dt), True)
secr = (Secr(prim = prim, dt = dt), False)
nucleo = (Nucleo(prim = prim, dt = dt), True)
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


gv[0].alim = False #LOHS

secr_on = False
scram = False
i_secr = -1
t_scram = -1
t_secr = -1

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
    
    if(new_P >= 13 and scram == False):
        nucleo[0].scram = True
        scram = True
        t_scram = prim.time[-1]
        s, f = nucleo
        s.t_scarm = t_scram
    if(new_P >= 13.7 and secr_on == False):
        prim.sis_on(secr)
        secr_on = True
        i_secr = i
        t_secr = prim.time[-1]
    
t0 = prim.time
p0 = prim.P




prim = primario()
estr = (Estruct(prim = prim, dt = dt), True)
secr = (Secr(prim = prim, dt = dt), False)
nucleo = (Nucleo(prim = prim, dt = dt, Q_th_nom=150E3), True)
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


gv[0].alim = False #LOHS

secr_on = False
scram = False
i_secr = -1
t_scram = -1
t_secr = -1

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
    
    if(new_P >= 13 and scram == False):
        nucleo[0].scram = True
        scram = True
        t_scram = prim.time[-1]
        s, f = nucleo
        s.t_scarm = t_scram
    if(new_P >= 13.7 and secr_on == False):
        prim.sis_on(secr)
        secr_on = True
        i_secr = i
        t_secr = prim.time[-1]

t1 = prim.time
p1 = prim.P








prim = primario()
estr = (Estruct(prim = prim, dt = dt), True)
secr = (Secr(prim = prim, dt = dt), False)
nucleo = (Nucleo(prim = prim, dt = dt, Q_th_nom=50E3), True)
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


gv[0].alim = False #LOHS

secr_on = False
scram = False
i_secr = -1
t_scram = -1
t_secr = -1

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
    
    if(new_P >= 13 and scram == False):
        nucleo[0].scram = True
        scram = True
        t_scram = prim.time[-1]
        s, f = nucleo
        s.t_scarm = t_scram
    if(new_P >= 13.7 and secr_on == False):
        prim.sis_on(secr)
        secr_on = True
        i_secr = i
        t_secr = prim.time[-1]

t2 = prim.time
p2 = prim.P




















fig3 = Figure(figsize=(12,12))
ax = fig3.add_subplot(1,1,1)
ax.plot(t0, p0, label="Q = 100 MW")
ax.plot(t1, p1, label="Q = 150 MW")
ax.plot(t2, p2, label="Q = 50 MW")
ax.grid(axis="both")
#ax.plot([t_scram, t_scram], [12,14], '--', label="SCRAM", color="black")
#ax.plot([t_secr, t_secr], [12,14], '--', label="inicio SECR", color="gray")
ax.set_title("Presion")
ax.set_ylabel("Presion [MPa]")
ax.set_xlabel("Tiempo [s]")
ax.legend()
fig3.show()




input()
#print(prim.Q_sis)