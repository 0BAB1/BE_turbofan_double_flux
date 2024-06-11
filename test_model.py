# Name : TEST FILE MODELE PARAMETRIQUE D'UN TURBOFAN DOUBLE FLUX
#
# BABIN-RIBY Hugo, BUFFARD-MICELI Eliot, PANZANI Mickael, ISSINA Fabert, BONNIN Baptiste
#
# File description : Test file to test the model with print and progressive reporting for dubug

from reporting import report
from math import sqrt
# données croisière + constantes
alt = 35e3 * 0.3048  # m
r = 287              # J.kg-1.K-1
rp = 291             # post combustion
gamma = 1.4
gammap = 1.33
P_k = 42800e3        # J.kg-1 pourvoir calorifique kero
LAMBDA = 11          # debit air / debit caburant
PI_f = 1.4           # taux de compresssion du fan

# données turbofan
XI_e = 0.98  # pertes de charges entrée d'air
ETA_c = 0.9  # rendement polytropique compresseur
ETA_f = 0.92  # rendement polytropique fan
ETA_comb = 0.99  # rendement de combustion
XI_cc = 0.95  # pertes de charges chambre de combustion
ETA_m = 0.98  # rendement arbre moteur
ETA_t_hp = 0.89  # rendement polytropique turbine HP
ETA_t_bp = 0.90  # rendement polytropique turbine BP
XI_tuy = 0.98  # perte de charges tuyere ejection
F_objectif = 21000  # N : force de poussée a atteindre

# temps & ratios
Tt4 = 1600  # K température totale de fin de combustion
OPR = 40  # pi_c = PT3/PT2
PI_c = OPR
PI_chp = 22  # taux de compression HP = Pt3/Pt2

# CONDITIONS EN 0
M0 = 0.78            # Mach
Ts0 = 217            # K
Ps0 = 22700          # Pa
Pt0 = Ps0 * (1 + (gamma-1)/2 * M0**2)**(gamma/(gamma-1))
Tt0 = Ts0 * (1 + (gamma-1) / 2 * M0**2)

report(0, Pt0, Tt0)

# CONDITIONS EN 2 => flux primaire

Pt2 = Pt0 * XI_e
Tt2 = Tt0
Ps2 = Ps0
Ts2 = Ts0

report(2, Pt2, Tt2)

# CONDITIONS EN 21 => Passage du fan

Pt21 = Pt2 * PI_f
Tt21 = Tt2 * PI_f ** ((gamma-1)/(gamma*ETA_f))

report(21, Pt21, Tt21)

# CONDITIONS EN 25 => Passage du compresseur BP

PI_cbp = OPR / (PI_f*PI_chp)
Pt25 = Pt21 * PI_cbp
Tt25 = Tt21 * PI_cbp ** ((gamma-1)/(gamma*1))

report(25, Pt25, Tt25)

# CONDITIONS EN 3 => Passage tous compresseurs

Pt3 = Pt25 * PI_chp
Tt3 = Tt25 * PI_chp ** ((gamma-1)/(gamma*ETA_c))

report(3, Pt3, Tt3)

# CONDITIONS EN 4 => passage chambre de combustion
Pt4 = Pt3 * XI_cc  # compression isobare
# Tt4 = 1600 K (données capteur)
report(4, Pt4, Tt4)

# CONDITIONS EN 5 => passage turbines, utilisation équilibre méca

# eq méca :    Pc = -Pt*ETA_m
#           => m'_c . Cp_c . DTt_c = -ETA_c . m'_t . Cp_t . DTt_t

# on determine les Cp & alpha
Cpc = (gamma*r)/(gamma-1)
Cpt = (gammap*rp)/(gammap-1)
alpha = (Tt4-Tt3)*Cpc/P_k

Tt45 = Tt4 + (Cpc * (Tt3-Tt25))/(-ETA_m*(1+alpha)*Cpt)
Pt45 = Pt4 * (Tt45/Tt4)**(gammap/((gammap-1)*ETA_t_hp))

report(45, Pt45, Tt45)

Tt5 = Tt45 + (Cpc * (Tt25-Tt21))/(-ETA_m*(1+alpha)*Cpt)
Pt5 = Pt45 * (Tt5/Tt45)**(gammap/((gammap-1)*ETA_t_bp))

report(5, Pt5, Tt5)

# CONDITIONS EN 9 => passage de la tuyère

Pt9 = Pt5 * XI_tuy
Tt9 = Tt5
Ps9 = Ps0
rho = Ps0/(r*Ts0)  # Kg/m3
v9 = sqrt((Pt9-Ps9*(2/rho)))
a = sqrt(gammap*rp*Ts0)
M9 = v9/a
print("En 9 => v9 : ", int(v9), " a : ", int(a), " Machs : ", (M9))

# CALCULS DES RENDEMENTS

v0 = a * M0
print("avec v0", v0)
Fp_spe = (1/(1+LAMBDA))*((1+alpha)*v9-v0)
print(Fp_spe*10, " N/(kg.s-1) fp spe")

# POUR LE FLUX SECONDAIRE

Fs_spe = (LAMBDA/(1+LAMBDA))*((1+alpha)*v9-v0)
print(Fs_spe*10, " N/(kg.s-1) fs spe")

# En total !

Ft_spe = Fs_spe + Fp_spe
debit_air = F_objectif/Ft_spe

print(debit_air, " kg.s-1")

# Calcul de la section d'entrée d'air

Rmoyeux = 0.3  # rapport moyeux rmin/rmax

Rmax = sqrt(debit_air/(rho*v0*3.14*(1-(Rmoyeux**2))))

print("on trouve rmax : ", Rmax,
      "donc D total moteur :", 2*Rmax*(1+Rmoyeux))

# Nb flux primaire validé : concordance avec data de Fabert

####################
# FLux secondiare
####################

# CONDITIONS EN 12
Pt12 = Pt0 * XI_e  # ok verifié aupres de data d'Eliot
Tt12 = Tt0
M12 = 0.6

# CONDITIONS EN 17
Pt17 = Pt12 * PI_f  # ok verifié aupres de data d'Eliot
Tt17 = Tt12 * (PI_f**((PI_f-1)/(PI_f*ETA_f)))

# conditions en 18 <=> conditions en 19 : on skip !
# CONDITIONS EN 19
Pt19 = Pt17 * XI_tuy  # ok verifié aupres de data d'Eliot
Tt19 = Tt17
Ps19 = Ps0
M19 = sqrt((((Pt19/Ps19)**((gamma-1)/gamma)) - 1)*2 /
           (gamma-1))  # ok verifié aupres de data d'Eliot
# v = M * a avec a = racine de gamma r t, ok verifié aupres de data d'Eliot
v19 = M19 * sqrt(gamma*r*Tt19)

# PUISSANCES ET RENDEMENTS
# PK : W/(kg/S)
# Ft_spe : N/(Kg/s)
debit_air_secondaire = LAMBDA * debit_air
P_cycle_primaire = 0.5*debit_air*(((1+alpha)*(v9**2))-(v0**2))
P_cycle_secondaire = 0.5*debit_air_secondaire*((v19**2)-(v0**2))
P_cycle = P_cycle_primaire + P_cycle_secondaire
P_chimique = debit_air*alpha*P_k
ETA_thermique = P_cycle/P_chimique
P_propulsive = F_objectif*v0
ETA_propulsif = P_propulsive / P_cycle
ETA_total = ETA_thermique * ETA_propulsif
print(ETA_thermique, ETA_propulsif, ETA_total)
