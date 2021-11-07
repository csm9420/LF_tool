# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:35:17 2020

@author: geuk_mi
"""

import CoolProp.CoolProp as CP

# CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'c:\\Program Files (x86)\\REFPROP\\')
# CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'RefpropEinbindung\\REFPROP\\')


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

