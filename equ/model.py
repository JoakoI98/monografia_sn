from .constants import *
import numpy as np

class primario:
    def __init__(self):
        #Definir variables del primario
        self.P = [prim_c.p_n]
        self.T = [T_sat(prim_c.p_n)]
        self.mg = [prim_c.vol_g*rhog(prim_c.p_n)]
        self.mf = [prim_c.vol_f*rhof(prim_c.p_n)]
        self.m = [self.mf[0] + self.mg[0]]
        self.VT = prim_c.vol_t
        self.vol_f = [prim_c.vol_f]
        self.vol_g = [prim_c.vol_g]
        self.time = [0]
        #Sistemas que intervienen
        self.sis = [] #Tupla con el sistema y un valor booleano que dice si esta prendido o no
        self.debug = True
        self.Q_sis = dict()
        self.Q_neta_tot = []
        
    
    @property
    def Q_neta(self):
        Q_n = 0
        for el, flag in self.sis:
            if flag:
                qi = el.Q
                if self.debug:
                    print(f"{type(el).__name__}     Q = {qi}")
                Q_n = Q_n + qi
                if not type(el).__name__ in self.Q_sis:
                    self.Q_sis[type(el).__name__] = []
                self.Q_sis[type(el).__name__].append(qi)
        self.Q_neta_tot.append(Q_n)
        return Q_n
        
    def print_sis(self):
        for el, flag in self.sis:
            if flag:
                print(f"ON {type(el).__name__}")
    
    def replace_sistem(self, new_sis):
        p = -1
        for i, el in enumerate(self.sis):
            if type(el[0]).__name__ == type(new_sis[0]).__name__:
                p = i
        if p != -1:
            self.sis[p] = new_sis
    
    def sis_on(self, new_sis):
        p = -1
        for i, el in enumerate(self.sis):
            if type(el[0]).__name__ == type(new_sis[0]).__name__:
                p = i
        if p != -1:
            s, f = self.sis[p]
            self.sis[p] = (s, True)
    
    def sis_off(self, new_sis):
        p = -1
        for i, el in enumerate(self.sis):
            if type(el[0]).__name__ == type(new_sis[0]).__name__:
                p = i
        if p != -1:
            s, f = self.sis[p]
            self.sis[p] = (s, False)
    
    @property
    def MP_neta(self):
        return sum([el.mp for el, f in self.sis if f == True])
    
    @property
    def t(self):
        return self.time[-1]
        
    def A(self, P):
        return (vg(P)*hf(P)-vf(P)*hg(P))/vfg(P)
    def B(self, P):
        return (hfg(P)/vfg(P))


class sys:
    def __init__(self, prim = None, dt = 0):
        self.h = self.A = self.c = self.m = 1
        self.prim = prim
        self.dt = dt
        self.debug = True
    
    @property
    def mp(self):
        return 0
    
    @property
    def Q(self):
        return None

class Estruct(sys):
    def __init__(self, prim = None, dt = 0):
        sys.__init__(self, prim, dt)
        self.h = 0.3 #kW/m2 K
        self.A = 300 #float(pi)*6*12*(mRPV1+mRPV2)/mRPV1 #m2
        self.T = [T_sat(P0)]
        self.m = (mRPV1+mRPV2)
        self.c = 0.51 #kJ/kg K
        
    @property
    def Q(self):
        Qi = self.h*self.A*(self.T[-1] - self.prim.T[-1])
        self.T.append(self.T[-1] - self.dt * Qi/(self.m * self.c))
        return Qi

class Secr(sys):
    def __init__(self, prim = None, dt = 0, q0 = 2E3):
        sys.__init__(self, prim, dt)
        self.T = [T_sat(Patm)]
        self.QSECR0 = q0 #kW
        self.h=self.QSECR0/(T_sat(P0)-self.T[0])
    
    @property
    def Q(self):
        Qi = self.h*self.A*(self.T[-1] - self.prim.T[-1])
        return Qi
    
class Nucleo(sys):
    def __init__(self, prim = None, dt = 0, Q_th_nom = 100E3):
        sys.__init__(self, prim, dt)
        Tcomb0=600 #grados Celsius
        self.Q0 = prim_c.pot
        r1=3.8e-3
        L=1.4
        Ncomb=6583
        v1=np.pi*r1**2*L*Ncomb
        r2=3.875e-3
        r3=4.5e-3
        v2=np.pi*(r3**2-r2**2)*L*Ncomb
        fv1=(v1)/(v1+v2)
        fv2=(v2)/(v1+v2)
        rhoc1=3e3
        rhoc2=2e3
        ccomb=(fv1*rhoc1+fv2*rhoc2)
        self.m=ccomb*(v1+v2)
        self.h = self.Q0/(Tcomb0-self.prim.T[0])
        self.T = [Tcomb0]
        self.scram = False
        self.Q_th_nom = Q_th_nom
        self.t_scarm = 0
    
    @property
    def Q(self):
        Qi = self.h*self.A*(self.T[-1] - self.prim.T[-1])
        if self.scram:
            delta_T = (q_dec(self.prim.t - self.t_scarm, self.Q0)-Qi)/(self.m*self.c)
        else:
            delta_T = (self.Q_th_nom-Qi)/(self.m*self.c)
        if self.debug:
            print(f"Nucleo:   Qi = {Qi} A = {self.A}   h = {self.h}    Q0 = {self.Q0} tcomb0 = {self.T[0]}  tp0 = {self.prim.T[0]}  T = {self.T[-1]}    Tp = {self.prim.T[-1]}    deltaT = {delta_T}")
        self.T.append(self.T[-1] + delta_T*self.dt)
        return Qi
    
class GV(sys):
    def __init__(self, prim = None, dt = 0):
        sys.__init__(self, prim, dt)
        PGV=4.7 #MPa
        TGV=T_sat(PGV)
        self.T = [TGV]
        QGV0=prim_c.pot
        self.h=QGV0/(self.prim.T[0]-TGV) #kW/K
        self.m0 = 600 #kg
        self.m = [self.m0]
        self.alim = True
        self.P = PGV
        self.Q0 = QGV0
    
    @property
    def Q(self):
        Qi = self.h*self.A*(self.m[-1]/self.m0)*(self.T[-1] - self.prim.T[-1])
        if self.alim == False:
            self.m.append(self.m[-1] + self.dt * Qi / hfg(self.P))
        else:
            self.m.append(self.m0)
        if self.debug:
            print(f"GV:   Qi = {Qi} A = {self.A}   h = {self.h}    Q0 = {self.Q0} tcomb0 = {self.T[0]}  tp0 = {self.prim.T[0]}  T = {self.T[-1]}    Tp = {self.prim.T[-1]}  m = {self.m[-1]}   m0 = {self.m0}")
        return Qi

class LOCA(sys):
    def __init__(self, prim = None, dt = 0):
        sys.__init__(self, prim, dt)
        dLOCA=1.5*2.54/100
        self.A=np.pi*(dLOCA**2)/4
        self.mp_t = []
    
    
    @property
    def mp(self):
        P = self.prim.P[-1]
        if P/Patm > 2:
            mp = self.A*np.sqrt(k*P*1e6*rhog(P)*((2/(k+1))**((k+1)/(k-1))))
        else:
            mp = self.A * np.sqrt(2*rhog(P)*(P-Patm)*1E6)
        return -1*mp
    
    @property
    def Q(self):
        self.mp_t.append(self.mp)
        return self.mp * hg(self.prim.P[-1])
    
class SIE(sys):
    def __init__(self, prim = None, dt = 0):
        sys.__init__(self, prim, dt)
        ASIE=4.53e-3
        kFSIE = 2.47e4

        PcritSIE=1.5 #Mpa
        hfSIE=188.44
        rhofSIE=1000 #kg/m3
        gamma=1.0

        RN2=0.2968
        VSIE0_N2 = 19.0 #m3, es el de N2
        PSIE0_N2 = 2.0 #MPa
        TSIE0_N2 = 45 #C
        mSIE_N2 = PSIE0_N2*VSIE0_N2/(RN2*(TSIE0_N2+273.15))

        VSIEf_N2 = 69.0
        PSIEf_N2 = PSIE0_N2*(VSIE0_N2/VSIEf_N2)**gamma
        TSIEf_N2 = (PSIEf_N2*VSIEf_N2/(RN2*mSIE_N2))-273.15
    
        self.P = [PSIE0_N2]
        self.vol_n2 = [VSIE0_N2]
        self.A = ASIE
        self.rho = rhofSIE
        self.k = kFSIE
        self.pmin = PSIEf_N2
        self.gamma = gamma
        self.hf = hfSIE
        self.mp_t = []
        
    @property
    def mp(self):
        P = self.prim.P[-1]
        mp = 0
        if self.P[-1] > self.pmin and self.P[-1] > P:
            mp = self.A*np.sqrt(2*self.rho*((self.P[-1]-P)*1e6)/self.k)
        return mp
    
    @property
    def Q(self):
        mp = self.mp
        self.vol_n2.append(self.vol_n2[-1] + mp/self.rho * self.dt)
        new_p = self.P[-1]*(self.vol_n2[-2]/self.vol_n2[-1])**self.gamma
        self.P.append(new_p)
        self.mp_t.append(mp)
        return mp * self.hf


class pot_cts(sys):
    def __init__(self, pot = 1, prim = None, dt = 0):
        sys.__init__(self, prim, dt)
        self.pot = pot
    
    @property
    def Q(self):
        return self.pot