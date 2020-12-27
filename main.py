from equ import *
from matplotlib.pyplot import figure as Figure


def d_dx(f,x0):
    dx=1e-4
    return (f(x0+dx)-f(x0-dx))/(2*dx)

N=50000# 48*3600
tTotal=36000# 48*3600.0
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
    
t = prim.time
estr, f= estr
nucleo, f= nucleo
gv, f= gv




fig2 = Figure(figsize=(14,12))
a0 = fig2.add_subplot(2,2,1)
a1 = fig2.add_subplot(2,2,2)
a2 = fig2.add_subplot(2,1,2)

a0.plot(t, prim.P, label="P Primario")
a1.plot(t, gv.m, label="Masa GV")
a2.plot(t, prim.T, label="T Primerio")
a2.plot(t, nucleo.T, label="T Nucleo")

a0.set_title("Presion")
a1.set_title("Masa")
a2.set_title("Temperatura")

a0.set_ylabel("Presion [MPa]")
a1.set_ylabel("Masa [Kg]")
a2.set_ylabel("Temperatura [°C]")

a0.set_xlabel("Tiempo [s]")
a1.set_xlabel("Tiempo [s]")
a2.set_xlabel("Tiempo [s]")


a0.grid(axis="both")
a1.grid(axis="both")
a2.grid(axis="both")

a0.plot([t_scram, t_scram], [12,14], '--', label="SCRAM", color="black")
a0.plot([t_secr, t_secr], [12,14], '--', label="inicio SECR", color="gray")

a1.plot([t_scram, t_scram], [0,620], '--', label="SCRAM", color="black")
a1.plot([t_secr, t_secr], [0,620], '--', label="inicio SECR", color="gray")

a2.plot([t_scram, t_scram], [310,600], '--', label="SCRAM", color="black")
a2.plot([t_secr, t_secr], [310,600], '--', label="inicio SECR", color="gray")

a0.legend()
a1.legend()
a2.legend()
fig2.show()



fig3 = Figure(figsize=(12,12))
ax = fig3.add_subplot(1,1,1)
ax.plot(t, prim.P, label="P Primario")
ax.grid(axis="both")
#ax.plot([t_scram, t_scram], [12,14], '--', label="SCRAM", color="black")
#ax.plot([t_secr, t_secr], [12,14], '--', label="inicio SECR", color="gray")
ax.set_title("Presion")
ax.set_ylabel("Presion [MPa]")
ax.set_xlabel("Tiempo [s]")
ax.legend()
fig3.show()

fig = Figure(figsize=(14,6))
axes = [fig.add_subplot(1,2,i) for i in range(1,3)]

serc_pot = np.concatenate([np.zeros(i_secr),np.asarray(prim.Q_sis['Secr'])])
t.remove(0)
axes[0].plot(t, np.asarray(prim.Q_sis['Nucleo'])*1E-3, label="Nucleo")
axes[0].plot(t, np.abs(serc_pot)*1E-3, label="SECR")
axes[0].plot(t, np.abs(np.asarray(prim.Q_sis['GV']))*1E-3, label="GV")
axes[0].plot([t_scram, t_scram], [0,1E3], 'b--', label="SCRAM", color="black")
axes[0].plot([t_secr, t_secr], [0,1E3], 'g--', label="inicio SECR", color="gray")


axes[1].plot(t, np.asarray(prim.Q_neta_tot)*1E-3, label="Q Neto")
axes[1].plot([t_scram, t_scram], [0,1E3], '--', label="SCRAM", color="black")
axes[1].plot([t_secr, t_secr], [0,1E3], '--', label="inicio SECR", color="gray")

for e in axes:
    e.grid(axis = 'both')
    e.legend()
    e.set_xlabel("tiempo [s]")
    e.set_ylabel("Potencia [MW]")
    e.set_yscale('log')
    e.set_ylim(1E-1, 1.1*1E2)

axes[0].set_title("Potencia")
axes[1].set_title("Potencia Neta")

fig.show()





input()
#print(prim.Q_sis)