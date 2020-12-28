from equ import *
from matplotlib.pyplot import figure as Figure


def d_dx(f,x0):
    dx=1e-4
    return (f(x0+dx)-f(x0-dx))/(2*dx)

N=20*3600# 48*3600
tTotal=50*3600# 48*3600.0
dt=tTotal/N




prim = primario()
estr = (Estruct(prim = prim, dt = dt), True)
secr = (Secr(prim = prim, dt = dt), False)
nucleo = (Nucleo(prim = prim, dt = dt), True)
gv = (GV(prim = prim, dt = dt), True)
loca = (LOCA(prim = prim, dt = dt), True)
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
nucleo[0].scram = True

secr_on = False
i_secr = -1
t_secr = -1

sie_on = False
i_sie = -1
t_sie = -1

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
    

    if((new_P >= 13.7 or new_P <= 10) and secr_on == False):
        prim.sis_on(secr)
        secr_on = True
        i_secr = i
        t_secr = prim.time[-1]
        print(f"SECR: {t_secr}")
    if( new_P <= 1.5 and sie_on == False):
        prim.sis_on(sie)
        sie_on = True
        i_sie = i
        t_sie = prim.time[-1]
        print(f"SIE: {t_sie}    {new_P}")
    
t = prim.time
estr, f= estr
nucleo, f= nucleo
gv, f= gv
loca, f= loca
sie, f= sie




#t.remove(0)
fig = Figure(figsize=(12,12))
axes0 = fig.add_subplot(1,1,1)
axes0.plot(np.asarray(t)/3600, prim.P, label="P primario (Inyeccion a 1.5MPa)")
axes0.plot([t_sie/3600, t_sie/3600], [0, 12.8],'--' ,label=f"Entrada del sie (Inyeccion a 1.5MPa) t={t_sie/3600:.2f} min", color="black")


m_sie = np.concatenate((np.zeros(i_sie+1), np.asarray(sie.mp_t)))

fig2 = Figure(figsize=(12,12))
axes2 = fig2.add_subplot(1,1,1)
axes2.plot(np.asarray(t)/3600, m_sie, label="M SIE (Inyeccion a 1.5MPa)")


fig3 = Figure(figsize=(12,12))
axes3 = fig3.add_subplot(1,1,1)
axes3.plot(np.asarray(t)/3600, prim.vol_f, label="Volumen liquido (Inyeccion a 1.5MPa)")












prim = primario()
estr = (Estruct(prim = prim, dt = dt), True)
secr = (Secr(prim = prim, dt = dt), False)
nucleo = (Nucleo(prim = prim, dt = dt), True)
gv = (GV(prim = prim, dt = dt), True)
loca = (LOCA(prim = prim, dt = dt), True)
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
nucleo[0].scram = True

secr_on = False
i_secr = -1
t_secr = -1

sie_on = False
i_sie = -1
t_sie = -1

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
    

    if((new_P >= 13.7 or new_P <= 10) and secr_on == False):
        prim.sis_on(secr)
        secr_on = True
        i_secr = i
        t_secr = prim.time[-1]
        print(f"SECR: {t_secr}")
    if( new_P <= 1.9 and sie_on == False):
        prim.sis_on(sie)
        sie_on = True
        i_sie = i
        t_sie = prim.time[-1]
        print(f"SIE: {t_sie}    {new_P}")
    
t = prim.time
estr, f= estr
nucleo, f= nucleo
gv, f= gv
loca, f= loca
sie, f= sie




#t.remove(0)
axes0.plot(np.asarray(t)/3600, prim.P, label="P primario (Inyeccion a 1.9MPa)")
axes0.plot([t_sie/3600, t_sie/3600], [0, 12.8],'--' ,label=f"Entrada del sie (Inyeccion a 1.9MPa): t={t_sie/3600:.2f} min", color="black")


m_sie = np.concatenate((np.zeros(i_sie+1), np.asarray(sie.mp_t)))

axes2.plot(np.asarray(t)/3600, m_sie, label="M SIE (Inyeccion a 1.9MPa)")



axes3.plot(np.asarray(t)/3600, prim.vol_f, label="Volumen liquido (Inyeccion a 1.9MPa)")



















prim = primario()
estr = (Estruct(prim = prim, dt = dt), True)
secr = (Secr(prim = prim, dt = dt), False)
nucleo = (Nucleo(prim = prim, dt = dt), True)
gv = (GV(prim = prim, dt = dt), True)
loca = (LOCA(prim = prim, dt = dt), True)
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
nucleo[0].scram = True

secr_on = False
i_secr = -1
t_secr = -1

sie_on = False
i_sie = -1
t_sie = -1

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
    

    if((new_P >= 13.7 or new_P <= 10) and secr_on == False):
        prim.sis_on(secr)
        secr_on = True
        i_secr = i
        t_secr = prim.time[-1]
        print(f"SECR: {t_secr}")
    if( new_P <= 1 and sie_on == False):
        prim.sis_on(sie)
        sie_on = True
        i_sie = i
        t_sie = prim.time[-1]
        print(f"SIE: {t_sie}    {new_P}")
    
t = prim.time
estr, f= estr
nucleo, f= nucleo
gv, f= gv
loca, f= loca
sie, f= sie




#t.remove(0)

axes0.plot(np.asarray(t)/3600, prim.P, label="P primario (Inyeccion a 1MPa)")
axes0.plot([t_sie/3600, t_sie/3600], [0, 12.8],'--' ,label=f"Entrada del sie (Inyeccion a 1MPa): t={t_sie/3600:.2f} min", color="black")



axes2.plot(np.asarray(t)/3600, np.concatenate((np.zeros(i_sie+1), np.asarray(sie.mp_t))), label="M SIE (Inyeccion a 1MPa)")


axes3.plot(np.asarray(t)/3600, prim.vol_f, label="Volumen liquido (Inyeccion a 1MPa)")





























#axes.plot(t, prim.Q_sis['GV'], label="Gv")
#axes.plot(t, prim.Q_sis['Estruct'], label="estr")
#axes.plot(t, prim.Q_sis['Nucleo'], label="core")
#axes.plot(t, prim.Q_sis['LOCA'], label="loca")
#axes.plot(t[i_secr:], prim.Q_sis['Secr'], label="secr")


#fig = Figure(figsize=(12,20))
#axes = [fig.add_subplot(2,2,i) for i in range(1,5)]
#axes[0].plot(t, prim.P)
#axes[1].plot(t, prim.T)
#axes[2].plot(t, estr.T)
#axes[3].plot(t, nucleo.T)
#
#axes[0].set_title("Presion en el primario")
#axes[1].set_title("Temperatura en el primario")
#axes[2].set_title("Temperatura en estructuras")
#axes[3].set_title("Temperatura en el combustible")
#
#for e in axes:
#    e.grid(axis = 'both')
#

    
#fig = Figure(figsize=(12,12))
#axes = fig.add_subplot(1,1,1) 
#t.remove(0)
#axes.plot(t, prim.Q_sis['GV'], label="Gv")
#axes.plot(t, prim.Q_sis['Estruct'], label="estr")
#axes.plot(t, prim.Q_sis['Nucleo'], label="core")
#axes.plot(t, prim.Q_sis['LOCA'], label="loca")
#axes.plot(t[i_secr:], prim.Q_sis['Secr'], label="secr")
#axes.legend()
#axes.grid(axis = 'both')












axes3.plot([t[0]/3600, t[-1]/3600], [prim_c.vol_nuc, prim_c.vol_nuc], label="Volumen de nucleo")


axes0.legend()
axes0.grid(axis = 'both')
axes0.set_xlabel("Tiempo [h]")
axes0.set_ylabel("Presion [MPa]")
axes0.set_title("Presion")

axes2.legend()
axes2.grid(axis = 'both')
axes2.set_xlabel("Tiempo [h]")
axes2.set_ylabel("Caudal [kg/s]")
axes2.set_title("SIE")

axes3.legend()
axes3.grid(axis = 'both')
axes3.set_xlabel("Tiempo [h]")
axes3.set_ylabel("Caudal [kg/s]")
axes3.set_title("Volumen de liquido")


fig.show()
fig2.show()
fig3.show()


input()
#print(prim.Q_sis)