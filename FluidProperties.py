# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 11:29:11 2020

@author: geuk_mi
"""

from Refprop import getRefPropValue
# from Refprop import getRefPropValue, getCoolPropValue
# import sys


def fluid_properties(fuel, T_iF, T_iO, p_cc):
    """ Calculating the fluid properties of the fuel and the oxidizer with RefProp

    :param fuel: Defining the type of fuel, 'H2' for hydrogen, 'CH4' for methane
    :type fuel: str
    :param T_iF: Temperature [K] of the fuel at the injector
    :type T_iF: float
    :param T_iO: Temperature [K] of the oxidizer at the injector
    :type T_iO: float
    :param p_cc: Theoretical pressure [bar] inside the combustion chamber
    :type p_cc: float

    :return: molar enthalpy of the fuel [J/kmol] and molar enthalpy of the oxidizer [J/kmol]
    :rtype : tuple[float, float]
    """

    # Note: As previous name of input parameter, p_cc_theo was defined, but p_cc was used. Because it makes sense to
    # use p_cc to calculate the enthalpies in the chamber, the name if the input parameter was changed accordingly.

    # TODO: The environment values are not really consistent for fuel and oxidizer
    # TODO: Why - 74873 [J/kmol] for fuel='CH4' ?

    # H_i = (H - H_env) * 1e-3 * M (with i = {O, F})
    # - H_i is enthalpy in J/kmol
    # - H is enthalpy in J/kg at T=T_iF [K] and p=p_cc_theo [bar]
    # - H_env is enthalpy in J/kg at environment conditions with T=298.15 [K] and p=101325 [Pa]
    # - M is molar mass in [kg/mol] at T=T_iF and p=100000 [Pa]

    # Oxidizer
    H_O = getRefPropValue('H', 'T', T_iO, 'P', p_cc * 100000, 'Oxygen') - getRefPropValue(
        'H', 'T', 298.15, 'P', 101325, 'Oxygen')
    H_O = H_O * getRefPropValue('M', 'T', 298.15, 'P', 101325, 'Oxygen') * 1e-3

    # Fuel
    # Check the fuel type first, 'H2' or 'CH4'
    if fuel == 'H2':
        H_F = getRefPropValue('H', 'T', T_iF, 'P', p_cc * 100000, 'Hydrogen') - getRefPropValue(
            'H', 'T', 298.15, 'P', 101325, 'Hydrogen')
        H_F = H_F * getRefPropValue('M', 'T', T_iF, 'P', 100000, 'Hydrogen') * 1e-3
    if fuel == 'CH4':
        H_F = getRefPropValue('H', 'T', T_iF, 'P', p_cc * 100000, 'Methane') - getRefPropValue(
            'H', 'T', 298.15, 'P', 101325, 'Methane')
        H_F = H_F * getRefPropValue('M', 'T', T_iF, 'P', 100000, 'Methane') * 1e-3 - 74873  # [kmol]

    # For a wrong fuel type, H_F will not be assigned but referenced in the return statement, this is intentional and
    # will crash the software if a unknown fuel type is given

    return H_F, H_O


# class FluidProperties:
#
#     def __init__(self, fuel, T_iF, T_iO, p_cc_theo):
#         self.fuel = fuel
#         self.TIF = T_iF
#         self.TIO = T_iO
#         self.p_cc_theo = p_cc_theo
#
#         # Fluid properties for fuel
#         if fuel == "H2":
#
#             # H_F = (H - H_env) * 1e-3 / M, whereby
#             # - H_F is enthalpy in J/kmol
#             # - H is enthalpy in J/kg at T=T_iF [K] and p=p_cc_theo [bar]
#             # - H_env is enthalpy in J/kg at environment conditions with T=298.15 [K] and p=101325 [Pa]
#             # - M is molar mass in [kg/mol] at T=T_iF and p=100000 [Pa]
#             self.H_F = getRefPropValue('H', 'T', T_iF, 'P', p_cc_theo * 100000, 'Hydrogen') - getRefPropValue(
#                 'H', 'T', 298.15, 'P', 101325, 'Hydrogen')
#             self.H_F = self.H_F * getRefPropValue('M', 'T', T_iF, 'P', 100000, 'Hydrogen') * 1e-3
#
#             # self.rho_F = getRefPropValue('D', 'T', T_iF, 'P', p_cc_theo * 100000, 'Hydrogen')
#             # self.c_F = getRefPropValue('A', 'T', T_iF, 'P', p_cc_theo * 100000, 'Hydrogen')
#
#         elif fuel == "CH4":
#             # TODO: Why - 74873 [J/kmol] ?
#             self.H_F = getRefPropValue('H', 'T', T_iF, 'P', p_cc_theo * 100000, 'Methane') - getRefPropValue(
#                 'H', 'T', 298.15, 'P', 101325, 'Methane')
#             self.H_F = self.H_F * getRefPropValue('M', 'T', T_iF, 'P', 100000, 'Methane') * 1e-3 - 74873
#
#             # self.rho_F = getRefPropValue('D', 'T', T_iF, 'P', p_cc_theo * 100000, 'Methane')
#             # self.c_F = getRefPropValue('A', 'T', T_iF, 'P', p_cc_theo * 100000, 'Methane')
#
#         else:
#             print("Error in enthalpy calculation. No valid fuel")
#
#         # Fluid properties for oxidizer
#         self.H_O = getRefPropValue('H', 'T', T_iO, 'P', p_cc_theo * 100000, 'Oxygen') - getRefPropValue(
#             'H', 'T', 298.15, 'P', 101325, 'Oxygen')
#         self.H_O *= getRefPropValue('M', 'T', 298.15, 'P', 101325, 'Oxygen') * 1e-3
#
#         # self.rho_O = getRefPropValue('D', 'T', T_iO, 'P', p_cc_theo * 100000, 'Oxygen')
#         # self.c_O = getRefPropValue('A', 'T', T_iO, 'P', p_cc_theo * 100000, 'Oxygen')
#
#     # def calculateFluidProperties(self):
#     #
#     #     if (self.fuel == "H2"):
#     #         self.H_F = getRefPropValue('H', 'T', self.TIF, 'P', self.p_cc_theo * 100000, 'Hydrogen') -
#     #           getRefPropValue('H', 'T', 298.15, 'P', 101325, 'Hydrogen')
#     #         self.H_F = self.H_F * getRefPropValue('M', 'T', self.TIF, 'P', 100000, 'Hydrogen') * 1e-3
#     #
#     #         self.rho_F = getRefPropValue('D', 'T', self.TIF, 'P', self.p_cc_theo * 100000, 'Hydrogen')
#     #         self.c_F = getRefPropValue('A', 'T', self.TIF, 'P', self.p_cc_theo * 100000, 'Hydrogen')
#     #
#     #     elif (self.fuel == "CH4"):
#     #         self.H_F = getRefPropValue('H', 'T', self.TIF, 'P', self.p_cc_theo * 100000, 'Methane') - getRefPropValue(
#     #             'H', 'T', 298.15, 'P', 101325, 'Methane')
#     #         self.H_F = self.H_F * getRefPropValue('M', 'T', self.TIF, 'P', 100000, 'Methane') * 1e-3 - 74873
#     #         self.rho_F = getRefPropValue('D', 'T', self.TIF, 'P', self.p_cc_theo * 100000, 'Methane')
#     #         self.c_F = getRefPropValue('A', 'T', self.TIF, 'P', self.p_cc_theo * 100000, 'Methane')
#     #     else:
#     #         print("Error in enthalpy calculation. No valid fuel")
#     #
#     #     self.H_O = getRefPropValue('H', 'T', self.TIO, 'P', self.p_cc_theo * 100000, 'Oxygen') - \
#     #                getRefPropValue('H', 'T', 298.15, 'P', 101325, 'Oxygen')
#     #     self.H_O *= getRefPropValue('M', 'T', 298.15, 'P', 101325, 'Oxygen') * 1e-3
#     #     self.rho_O = getRefPropValue('D', 'T', self.TIO, 'P', self.p_cc_theo * 100000, 'Oxygen')
#     #     self.c_O = getRefPropValue('A', 'T', self.TIO, 'P', self.p_cc_theo * 100000, 'Oxygen')
