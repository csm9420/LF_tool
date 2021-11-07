# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 11:07:09 2020

@author: geuk_mi
"""

import math
import sys
from FluidProperties import fluid_properties
import CEAFunctions
from Refprop import getRefPropValue, getCoolPropValue


class CombustionChamber:

    def __init__(self, config, d_cc, l_cc, epsilon_cc, p_cc, eta_cc, rof, mDot_O2, mDot_fuel_primary, fuel, T_if, T_io):
        """ Creates a combustion chamber.

        :param config:
        """

        # Geometry
        self.diameter = d_cc * 1e-3  # Diameter [m]
        self.length = l_cc * 1e-3  # Length [m]
        self.epsilon_cc = epsilon_cc  # Contraction ratio [-]
        self.area = math.pi / 4 * math.pow(self.diameter, 2)  # Cross section area of the chamber [m^2]

        # Load point definition
        self.p = p_cc  # Pressure [bar]
        self.eta_cc = eta_cc  # Combustion efficiency [-]
        self.rof = rof  # Rate of fuel (?) at primary injection [-]
        self.mDot_O2 = mDot_O2  # Mass flow oxidizer [kg/s]
        self.mDot_F_primary_injector = mDot_fuel_primary  # Primary fuel mass flow [kg/s]
        self.p_theo = self.p / self.eta_cc  # Theoretical pressure (nozzle, mass flow) [bar]

        # Fluid properties
        # H_F is fuel enthalpy [J/kmol]
        # H_O is oxidizer enthalpy [J/kmol]
        self.H_F, self.H_O = fluid_properties(fuel=fuel, T_iF=T_if, T_iO=T_io, p_cc=self.p)

        # Combustion temperatures
        self.T_theo = CEAFunctions.ceaIdealCombustionTemp(fuel=fuel, T_iF=T_if, h_F=self.H_F, T_iO=T_io,
                                                          h_O=self.H_O, rof=self.rof, p_cc=self.p,
                                                          accat=self.epsilon_cc)  # Theoretical temperature [K]
        self.T_combustion = self.T_theo * math.pow(self.eta_cc, 2)  # Combustion temperature [K]

        # Density
        # ceaRealCharacteristics returns a vector of values, but only the first element [0] is needed
        self.rho = CEAFunctions.ceaRealCharacteristics(fuel=fuel, T_iF=T_if, h_F=self.H_F, T_iO=T_io, h_O=self.H_O,
                                                       T_cc=self.T_combustion, rof=self.rof, p_cc=self.p,
                                                       accat=self.epsilon_cc)[0]

        # MWg (?)
        self.MWg_Chamber = CEAFunctions.ceaGetMoleMassProduct(fuel=fuel, T_iF=T_if, h_F=self.H_F, T_iO=T_io,
                                                              h_O=self.H_O, rof=self.rof, p_cc=self.p,
                                                              accat=self.epsilon_cc)
        self.MWg_Throat = CEAFunctions.ceaGetMoleMassProduct(fuel=fuel, T_iF=T_if, h_F=self.H_F, T_iO=T_io,
                                                             h_O=self.H_O, rof=self.rof, p_cc=self.p,
                                                             accat=self.epsilon_cc, throat=True)

        # TODO: replace for better usage, if time, this is not urgent
        self.config = config
        self.calculateChamberParameters()

    def calculateChamberParameters(self):

        # self.diameter = self.config.d_cc * 1e-3
        # self.length = self.config.l_cc * 1e-3
        # self.contractionRatio = self.config.epsilon_c
        # self.area = math.pi / 4 * math.pow(self.diameter, 2)

        # self.p = self.config.p_cc
        # self.combustionEfficiency = self.config.eta_cc
        # self.p_theo = self.p / self.combustion_efficiency
        # self.rof = self.config.rof_cc
        # self.mDot_O2 = self.config.mDot_O2
        # self.mDot_Fuel_PrimInj = self.config.mDot_F

        # fluidProp = fluidProperties(self.config.fuel, self.config.T_if, self.config.T_io, self.p)
        # self.fuelEnthalpie = fluidProp.H_F
        # self.oxidizerEnthalpie = fluidProp.H_O

        # self.T_theo = CEAFunctions.ceaIdealCombustionTemp(self.config.fuel, self.config.T_if, self.fuelEnthalpie,
        #                                                   self.config.T_io, self.oxidizerEnthalpie, self.rof, self.p,
        #                                                   self.contraction_ratio)
        # self.T_combustion = self.T_theo * math.pow(self.combustion_efficiency, 2)
        # solutionVec = CEAFunctions.ceaRealCharacteristics(self.config.fuel, self.config.T_if, self.H_F,
        #                                                   self.config.T_io, self.H_O, self.T_combustion,
        #                                                   self.rof, self.p, self.epsilon_cc)
        # self.rho = solutionVec[0]

        # self.MWg_Chamber = CEAFunctions.ceaGetMoleMassProduct(self.config.fuel, self.config.T_if, self.H_F,
        #                                                       self.config.T_io, self.H_O, self.rof,
        #                                                       self.p, self.epsilon_cc)
        # self.MWg_Throat = CEAFunctions.ceaGetMoleMassProduct(self.config.fuel, self.config.T_if, self.H_F,
        #                                                      self.config.T_io, self.H_O, self.rof, self.p,
        #                                                      self.epsilon_cc, True)

        #        self.meanVelocity = (self.mDot_O2 + self.mDot_O2/self.rof) / (self.rho * math.pi / 4 * self.area)
        #        self.l_star = self.area * self.length / self.contractionRatio
        #        self.theta = self.length / self.meanVelocity

        if ((not math.isnan(self.config.mDot_SF)) or (not math.isnan(self.config.mDot_WC))):
            if ((not math.isnan(self.config.mDot_SF)) and math.isnan(self.config.mDot_WC)):
                self.T_additinalInjection = self.config.T_if
                self.mDot_fuel_additional = self.config.mDot_SF
            elif ((not math.isnan(self.config.mDot_WC)) and math.isnan(self.config.mDot_SF)):
                self.T_additinalInjection = self.config.T_wc
                self.mDot_fuel_additional = self.config.mDot_WC
            else:
                if (self.config.fuel == "CH4"):
                    fld = "Methane"
                elif (self.config.fuel == "H2"):
                    fld = "Hydrogen"

                self.h_mix = (self.config.mDot_SF * getRefPropValue('H', 'T', self.config.T_if, 'P', self.p * 100000, fld)
                              + self.config.mDot_WC * getRefPropValue('H', 'T', self.config.T_wc, 'P',  self.p * 100000, fld)) / (self.config.mDot_SF + self.config.mDot_WC)

                self.T_additinalInjection = getRefPropValue('T', 'H', self.h_mix, 'P', self.p * 100000, fld)
                self.mDot_fuel_additional = self.config.mDot_SF + self.config.mDot_WC
        else:
            self.mDot_fuel_additional = 0
            self.T_additinalInjection = 0

        self.mDot_fuel = self.mDot_F_primary_injector + self.mDot_fuel_additional
        self.mDot = self.mDot_O2 + self.mDot_fuel

        self.drof = 0.1

        # Reading Combustion Chamber temperatures
        self.T_theo_plus = CEAFunctions.ceaIdealCombustionTempExtended(self.config.fuel, self.config.T_if,
                                                                       self.H_F, self.config.T_io,
                                                                       self.H_O, self.rof + self.drof,
                                                                       self.p, self.epsilon_cc)
        self.T_combustion_plus = [i * math.pow(self.eta_cc, 2) for i in self.T_theo_plus]
        self.T_theo_minus = CEAFunctions.ceaIdealCombustionTempExtended(self.config.fuel, self.config.T_if,
                                                                        self.H_F, self.config.T_io,
                                                                        self.H_O, self.rof - self.drof,
                                                                        self.p, self.epsilon_cc)
        self.T_combustion_minus = [i * math.pow(self.eta_cc, 2) for i in self.T_theo_minus]

        self.dt_drof_nozzleEntrance = (self.T_combustion_plus[1] - self.T_combustion_minus[1]) / (2 * self.drof)
        self.dt_drof_combustionFront = (self.T_combustion_plus[0] - self.T_combustion_minus[0]) / (2 * self.drof)

        # Reading c* and dc*dROF
        self.cstar = CEAFunctions.ceaGetCStar(self.config.fuel, self.config.T_if, self.H_F, self.config.T_io,
                                              self.H_O, self.rof, self.p, self.epsilon_cc)
        self.cstar_plus = CEAFunctions.ceaGetCStar(self.config.fuel, self.config.T_if, self.H_F,
                                                   self.config.T_io, self.H_O, self.rof + self.drof,
                                                   self.p, self.epsilon_cc)
        self.cstar_minus = CEAFunctions.ceaGetCStar(self.config.fuel, self.config.T_if, self.H_F,
                                                    self.config.T_io, self.H_O, self.rof - self.drof,
                                                    self.p, self.epsilon_cc)
        self.dc_drof = (self.cstar_plus - self.cstar_minus) / (2 * self.drof)

        # Reading chamber sound speed
        self.soundSpeed, self.machNumber, self.meanVelocity = CEAFunctions.ceaGetChamberSpeedValues(self.config.fuel,
                                                                                                    self.config.T_if,
                                                                                                    self.H_F,
                                                                                                    self.config.T_io,
                                                                                                    self.H_O,
                                                                                                    self.rof, self.p,
                                                                                                    self.epsilon_cc)
        self.l_star = self.area * self.length / self.epsilon_cc
        self.theta = self.length / self.meanVelocity

        # Reading specific heat ratio
        self.gamma = CEAFunctions.ceaGetSpecificHeatRatio(self.config.fuel, self.config.T_if, self.H_F,
                                                          self.config.T_io, self.H_O, self.rof, self.p,
                                                          self.epsilon_cc)
