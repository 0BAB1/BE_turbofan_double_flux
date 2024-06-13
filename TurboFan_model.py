# Name : MODELE PARAMETRIQUE D'UN TURBOFAN DOUBLE FLUX
#
# BABIN-RIBY Hugo, BUFFARD-MICELI Eliot, PANZANI Mickael, ISSINA Fabert
#
# File description : to plot data for turbo fan engine

# dont't need reporting here
# from reporting import report
from math import sqrt


class TurboFan():
    def __init__(self, OPR, Tt4, BPR=11, PIf=1.4) -> None:
        # données croisière + constantes
        self.alt = 35e3 * 0.3048  # m
        self.r = 287              # J.kg-1.K-1
        self.rp = 291             # post combustion
        self.gamma = 1.4
        self.gammap = 1.33
        self.P_k = 42800e3        # J.kg-1 pourvoir calorifique kero
        self.LAMBDA = BPR          # debit air / debit caburant
        self.PI_f = PIf           # taux de compresssion du fan
        # données turbofan
        self.XI_e = 0.98  # pertes de charges entrée d'air
        self.ETA_c = 0.9  # rendement polytropique compresseur
        self.ETA_f = 0.92  # rendement polytropique fan
        self.ETA_comb = 0.99  # rendement de combustion
        self.XI_cc = 0.95  # pertes de charges chambre de combustion
        self.ETA_m = 0.98  # rendement arbre moteur
        self.ETA_t_hp = 0.89  # rendement polytropique turbine HP
        self.ETA_t_bp = 0.90  # rendement polytropique turbine BP
        self.XI_tuy = 0.98  # perte de charges tuyere ejection
        # temps & ratios
        self.Tt4 = Tt4  # K température totale de fin de combustion
        self.OPR = OPR  # pi_c = PT3/PT2
        self.PI_c = self.OPR
        self.PI_chp = 22  # taux de compression HP = Pt3/Pt2
        self.F_objectif = 21000

        # CONDITIONS EN 0
        self.M0 = 0.78            # Mach
        self.Ts0 = 217            # K
        self.Ps0 = 22700          # Pa
        self.Pt0 = self.Ps0 * (1 + (self.gamma-1)/2 *
                               self.M0**2)**(self.gamma/(self.gamma-1))
        self.Tt0 = self.Ts0 * (1 + (self.gamma-1) / 2 * self.M0**2)

        # CONDITIONS EN 2 => flux primaire
        self.Pt2 = self.Pt0 * self.XI_e
        self.Tt2 = self.Tt0
        self.Ps2 = self.Ps0
        self.Ts2 = self.Ts0

        # CONDITIONS EN 21 => Passage du fan
        self.Pt21 = self.Pt2 * self.PI_f
        self.Tt21 = self.Tt2 * \
            self.PI_f ** ((self.gamma-1)/(self.gamma*self.ETA_f))

        # CONDITIONS EN 25 => Passage du compresseur BP
        self.PI_cbp = self.OPR / (self.PI_f*self.PI_chp)
        self.Pt25 = self.Pt21 * self.PI_cbp
        self.Tt25 = self.Tt21 * self.PI_cbp ** ((self.gamma-1)/(self.gamma*1))

        # CONDITIONS EN 3 => Passage tous compresseurs
        self.Pt3 = self.Pt25 * self.PI_chp
        self.Tt3 = self.Tt25 * \
            self.PI_chp ** ((self.gamma-1)/(self.gamma*self.ETA_c))

        # CONDITIONS EN 4 => passage chambre de combustion
        self.Pt4 = self.Pt3 * self.XI_cc  # compression isobare
        # Tt4 = 1600 K (données capteur)

        # CONDITIONS EN 5 => passage turbines, utilisation équilibre méca
        # eq méca :    Pc = -Pt*ETA_m
        #           => m'_c . Cp_c . DTt_c = -ETA_c . m'_t . Cp_t . DTt_t
        # on determine les Cp & alpha
        self.Cpc = (self.gamma*self.r)/(self.gamma-1)
        self.Cpt = (self.gammap*self.rp)/(self.gammap-1)
        self.alpha = 0.02  # This fixes alpha (more realistic)
        # self.Cpt = (self.gammap*self.rp)/(self.gammap-1) # This deduce alpha from Tt4
        self.Tt45 = self.Tt4 + (self.Cpc * (self.Tt3-self.Tt25)) / \
            (-self.ETA_m*(1+self.alpha)*self.Cpt)
        self.Pt45 = self.Pt4 * \
            (self.Tt45/self.Tt4)**(self.gammap/((self.gammap-1)*self.ETA_t_hp))

        self.Tt5 = self.Tt45 + \
            (self.Cpc * (self.Tt25-self.Tt21)) / \
            (-self.ETA_m*(1+self.alpha)*self.Cpt)
        self.Pt5 = self.Pt45 * \
            (self.Tt5/self.Tt45)**(self.gammap/((self.gammap-1)*self.ETA_t_bp))

        # CONDITIONS EN 9 => passage de la tuyère

        self.Pt9 = self.Pt5 * self.XI_tuy
        self.Tt9 = self.Tt5
        self.Ps9 = self.Ps0
        self.rho = self.Ps0/(self.r*self.Ts0)  # Kg/m3
        self.v9 = sqrt((self.Pt9-self.Ps9*(2/self.rho)))
        self.a = sqrt(self.gammap*self.rp*self.Ts0)
        self.M9 = self.v9/self.a

        ####################
        # FLux secondiare
        ####################

        # CONDITIONS EN 12
        self.Pt12 = self.Pt0 * self.XI_e  # ok vérifié auprès de data d'Eliot
        self.Tt12 = self.Tt0
        self.M12 = 0.6

        # CONDITIONS EN 17
        self.Pt17 = self.Pt12 * self.PI_f  # ok vérifié auprès de data d'Eliot
        self.Tt17 = self.Tt12 * \
            (self.PI_f**((self.PI_f-1)/(self.PI_f*self.ETA_f)))

        # conditions en 18 <=> conditions en 19 : on skip !
        # CONDITIONS EN 19
        self.Pt19 = self.Pt17 * self.XI_tuy  # ok vérifié auprès de data d'Eliot
        self.Tt19 = self.Tt17
        self.Ps19 = self.Ps0
        self.M19 = sqrt((((self.Pt19/self.Ps19)**((self.gamma-1)/self.gamma)) - 1)*2 /
                        (self.gamma-1))  # ok vérifié auprès de data d'Eliot
        # v = M * a avec a = racine de gamma r t, ok vérifié auprès de data d'Eliot
        self.v19 = self.M19 * sqrt(self.gamma*self.r*self.Tt19)

        # CALCULS DES POUSSEES
        self.v0 = self.a * self.M0
        self.Fp_spe = (1/(1+self.LAMBDA))*((1+self.alpha)*self.v9-self.v0)

        # POUR LE FLUX SECONDAIRE
        self.Fs_spe = (self.LAMBDA/(1+self.LAMBDA)) * \
            ((1+self.alpha)*self.v19-self.v0)

        # En total !
        self.Ft_spe = self.Fs_spe + self.Fp_spe
        self.debit_air = 21000 / \
            ((1+self.alpha)*(self.v9+self.LAMBDA*self.v19)-(self.LAMBDA+1)*self.v0)

        # Calcul de la section d'entrée d'air
        self.Rmoyeux = 0.3  # rapport moyeux rmin/rmax
        self.Rmax = sqrt(
            self.debit_air/(self.rho*self.v0*3.14*(1-(self.Rmoyeux**2))))

        # PUISSANCES ET RENDEMENTS
        # PK : W/(kg/S)
        # Ft_spe : N/(Kg/s)
        self.debit_air_secondaire = self.LAMBDA * self.debit_air
        self.P_cycle_primaire = 0.5*self.debit_air * \
            (((1+self.alpha)*(self.v9**2))-(self.v0**2))
        self.P_cycle_secondaire = 0.5 * \
            self.debit_air_secondaire*((self.v19**2)-(self.v0**2))

        self.P_cycle = self.P_cycle_primaire + self.P_cycle_secondaire
        self.P_chimique = self.debit_air*self.alpha*self.P_k
        self.ETA_thermique = self.P_cycle/self.P_chimique

        self.P_propulsive = self.F_objectif*self.v0
        self.ETA_propulsif = self.P_propulsive / self.P_cycle
        self.ETA_total = self.ETA_thermique * self.ETA_propulsif


engine = TurboFan(40, 1600, 11)
print(engine.ETA_thermique)
print(engine.Ft_spe, engine.Fp_spe, engine.Fs_spe)
