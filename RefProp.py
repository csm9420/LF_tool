# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:35:17 2020

@author: geuk_mi
"""

import CoolProp.CoolProp as CP

# CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'c:\\Program Files (x86)\\REFPROP\\')
CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'RefpropEinbindung\\REFPROP\\')


# Method handles communication with REFPROP and gives back the requested value
#   Current values:
#       reqVal: Requested value type:
#           A:  Speed of sound          [m/s]
#           C:  Cp                      [J/(kg K)]
#           D:  Density                 [kg/m^3]
#           H:  Enthalpy                [J/kg]
#           P:  Pressure                [Pa]
#           T:  Temperature             [K]
#           V:  Dynamic Viscosity       [Pa s]
#           I:  Surface Tension         [N/m]
#
#       spez1:  first input type:   T, P, H, D, C, R, M
#           T, P, H, D: see above
#           C:  Properties at the critical point
#           R:  Properties at the triple point
#           M:  properties at T_max and P_max
#
#           See http://www.coolprop.org/coolprop/HighLevelAPI.html#table-of-string-inputs-to-propssi-function


def getRefPropValue(reqVal, spez1, value1, spez2, value2, fluid):
    RefFluid = fluid  # "REFPROP::" + fluid
    return CP.PropsSI(reqVal, spez1, value1, spez2, value2, RefFluid)


def getCoolPropValue(reqVal, spez1, value1, spez2, value2, fluid):
    return CP.PropsSI(reqVal, spez1, value1, spez2, value2, fluid)


# This can be used for debug purposes
#
if __name__ == "__main__":

    # p = getRefPropValue('P', 'T', 99.914, 'D', 1049.7846, 'REFPROP::Oxygen')
    cp_ref = getRefPropValue('Cpmass', 'D', 42.7348, 'P', 565070.912, 'REFPROP::Hydrogen')
    cv_ref = getRefPropValue('Cvmass', 'D', 42.7348, 'P', 565070.912, 'REFPROP::Hydrogen')
    cp_cool = getRefPropValue('Cpmass', 'D', 42.7348, 'P', 565070.912, 'Hydrogen')
    cv_cool = getRefPropValue('Cvmass', 'D', 42.7348, 'P', 565070.912, 'Hydrogen')
    print(cp_ref)
    print(cp_cool)
    print('\n')
    print(cp_ref / cv_ref)



    # import matplotlib.pyplot as plt

    # rho = 1036
    # kappa = []
    # p_all = []
    # for p in range(int(1 * 1e5), int(20*1e5), 1000):
    #     p_all.append(p)
    #     kappa.append(getRefPropValue('ISENTROPIC_EXPANSION_COEFFICIENT', 'D', rho, 'P', p, 'Oxygen'))
    #
    # plt.figure()
    # plt.plot([x * 1e-5 for x in p_all], kappa)

    # T = []
    # rho =[]
    # for t in range(50, 500, 10):
    #     # rho.append(getCoolPropValue("V", "T", t, "P", 65 * 100000, "Oxygen") / getCoolPropValue("D", "T", t, "P", 65 * 100000, "Oxygen"))
    #     rho.append(getCoolPropValue("D", "T", t, "P", 65 * 100000, "Oxygen"))
    #     T.append(t)
    #
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(T, rho)
    # plt.xlabel('T [K]')
    # plt.ylabel('rho [kg\m^3]')
    # plt.title('Oxygen density, p = 65\u2009bar')
    # plt.grid()
    # plt.show()

    # from sympy import symbols, solve, sqrt, log
    # from scipy.optimize import newton_krylov
    #
    # Re = 10 ^ 6
    # d = 0.004
    # # epsilon = symbols('epsilon')
    # kw = 0.5
    # f = lambda epsilon: 1/sqrt(epsilon) + 2 * log(2.51 / (sqrt(epsilon))) - 0.27 * d / kw
    # sol = newton_krylov(lambda epsilon: 1/sqrt(epsilon) + 2 * log(2.51 / (sqrt(epsilon) * Re), 10) - 0.27 * d / kw, 0.25)
    # epsilon = lambda d, Re, kw: (float(newton_krylov(lambda epsilon: 1/sqrt(epsilon) + 2 * log(2.51 / (sqrt(epsilon) * Re), 10) - 0.27 * d / kw, 0.25)))

    # print(epsilon(d, Re, kw))

#
#     # For debug purposes, get some values
#     H_CH4_env = getRefPropValue('H', 'T', 298.15, 'P', 101325, 'Methane')  # J/kg
#     M_CH4 = getRefPropValue('M', 'T', 298.15, 'P', 100000, 'Methane') * 10e-3  # kg/kmol
#     h_CH4 = H_CH4_env * M_CH4
#
#     print('H_CHV_env = {}'.format(H_CH4_env))
#     print('M_CH4 = {}'.format(M_CH4))
#     print('h_CH4 = {}'.format(h_CH4))
