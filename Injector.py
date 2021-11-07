# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 08:48:46 2020

@author: geuk_mi
"""
import sys
import math
import Functions as func
from Refprop import getRefPropValue, getCoolPropValue
import numpy as np
import pandas as pd
import os

# Defining the directory where the injector configuration files are stored
FILE_DIRECTORY = 'Geometry_Injection_Configurations'

# Defining the ranges (= the rows) in the excel injector configuration file
# Do not change this, unless the excel file is changed
RANGE_LOX = range(1, 7)
RANGE_PRIMARY_FUEL = range(8, 14)
RANGE_SECONDARY_FUEL = range(15, 21)
RANGE_WINDOW_COOLING = range(22, 28)


class Element:

    def __init__(self, data, T, mDot):
        """ Represents a single manifold, orifice or line element.

        :param data: A :class:pandas.Series object containing all the information about a single element.
        """

        self.type = data.loc['type']
        self.length = data.loc['length']
        self.diameter = data.loc['hydraulic diameter']
        self.area = data.loc['area']
        self.volume = data.loc['volume']
        self.numOfElements = data.loc['num_elem']
        self.T = T
        self.mDot = mDot / self.numOfElements

        self.dp = 0  # pressure loss over element [bar]


class LOX:

    def __init__(self, data, T, mDot):
        """ Represents the LOX supply withing the injector.

        :param data : A
        """
        self.num_elements = len(data.columns)
        self.elements = [Element(data[k], T, mDot) for k in data]
        self.type = 'LOX'


class Fuel:
    def __init__(self, data, fluid_type, T, mDot):
        """ Represents the fuel supply withing the injector.

        :param data : A
        """
        self.num_elements = len(data.columns)
        self.elements = [Element(data[k], T, mDot) for k in data]
        self.type = fluid_type


class SecFuel:
    def __init__(self, data, fluid_type, T, mDot):
        """ Represents the secondary fuel supply withing the injector.

        :param data : A
        """
        self.num_elements = len(data.columns)
        self.elements = [Element(data[k], T, mDot) for k in data]
        self.type = fluid_type


class WinCool:
    def __init__(self, data, fluid_type, T, mDot):
        """ Represents the window cooling withing the injector/chamber.

        :param data : A
        """
        self.num_elements = len(data.columns)
        self.elements = [Element(data[k], T, mDot) for k in data]
        self.type = fluid_type


class injector_configuration:

    def __init__(self, config, injector_filename, fluid_fuel, fluid_window_cooling, T_io, T_if, T_sf, T_wc, mDot_O2,
                 mDot_F, mDot_SF, mDot_WC):
        """ Creates an injector.

        :param injector_file:
        :param config:
        """

        # Attributes
        self._filename = injector_filename

        # Open the excel injector file as pd.DataFrame;
        # - then split the DataFrame into four smaller pd.DataFrames to represent the four input areas.
        # The ranges in the excel file (= the rows) are hard-coded (see above) because the should never change;
        # - then rename the rows (index) inside the DataFrame based on '# Identifier' to have better access to the data
        # - then drop/remove the column '# Identifier' because it is not longer needed
        # - then drop/remove all empty columns (= columns containing only nan values) because they're not useful
        # The result are four DataFrames only containing the data, nothing else
        data = pd.read_excel(os.path.join(os.getcwd(), FILE_DIRECTORY, self._filename + '.xlsx'), engine='openpyxl')
        data_LOX = data.iloc[RANGE_LOX].rename(index=lambda k: data.loc[k]['# Identifier']) \
            .drop('# Identifier', axis='columns').dropna(axis='columns', how='all')
        data_primary_fuel = data.iloc[RANGE_PRIMARY_FUEL].rename(index=lambda k: data.loc[k]['# Identifier']) \
            .drop('# Identifier', axis='columns').dropna(axis='columns', how='all')
        data_secondary_fuel = data.iloc[RANGE_SECONDARY_FUEL].rename(index=lambda k: data.loc[k]['# Identifier']) \
            .drop('# Identifier', axis='columns').dropna(axis='columns', how='all')
        data_window_cooling = data.iloc[RANGE_WINDOW_COOLING].rename(index=lambda k: data.loc[k]['# Identifier']) \
            .drop('# Identifier', axis='columns').dropna(axis='columns', how='all')

        if not data_LOX.empty:
            self.LOX = LOX(data_LOX, T_io, mDot_O2)
        if not data_primary_fuel.empty:
            self.Fuel = Fuel(data_primary_fuel, fluid_fuel, T_if, mDot_F)
        if not data_secondary_fuel.empty:
            self.SecFuel = SecFuel(data_secondary_fuel, fluid_fuel, T_sf, mDot_SF)
        if not data_window_cooling.empty:
            self.WinCool = WinCool(data_window_cooling, fluid_window_cooling, T_wc, mDot_WC)

        # Previous code
        # self.injector_file = injector_filename
        self.config = config

        # self.parse_injector_config()  # Done
        # self.assignFluidTypes()  # Done
        # self.assignTemperaturesAndMassflowsToElements()  # Done

        self.recalculate_injector()

        # self.calculatePressureInElements()
        # self.calculateAcousticParameters()
        # if self.config.casMode:
        #     self.reduceInjectorHeadToSingleInjector()

    def recalculate_injector(self):
        """ Calculates injector parameters and assigned them to the injector.

        :rtype: None
        """
        # if self.config.use_whole_feed_system is True:
        #    # Includes the whole feed system
        #     self.calculate_feed_systems_parameter()
        #     injector_plots(os.path.join(self.config.result_path, 'Injector'), 'Oxygen_calculated.csv',
        #                   'Oxygen_scaled.csv', 'Oxygen_pressure.png', self.config.p_ot, self.config.p_io,
        #                   self.config.p_cc)
        #     injector_plots(os.path.join(self.config.result_path, 'Injector'), 'Hydrogen_calculated.csv',
        #                   'Hydrogen_scaled.csv', 'Hydrogen_pressure.png', self.config.p_ft, self.config.p_if,
        #                   self.config.p_cc)
        # else:
        #    # Approach by Geuking
        self.calculatePressureInElements()
        self.calculateAcousticParameters()
        if self.config.casMode is True:
           self.reduceInjectorHeadToSingleInjector()

        # Write parameter to file
        # with open('tools_for_post_processing/dp/injector_old2.csv', 'w') as file:
        #     import csv
        #     writer = csv.DictWriter(file, fieldnames=['type', 'p', 'dp', 'mdot', 'rho', 'v', 'T'])
        #     writer.writeheader()
        #     for element in self.LOX.elements:
        #         try:
        #             writer.writerow({'type': element.type, 'p': element.p, 'dp': element.dp, 'mdot': element.mDot,
        #                          'rho': element.rho, 'T': element.T, 'v': element.velocity,})
        #         except:
        #             # file.write('upps\n')
        #             writer.writerow({'type': element.type, 'p': element.p, 'dp': element.dp, 'mdot': element.mDot,
        #                              'rho': element.rho, 'T': element.T})

    # def parse_injector_config(self):
    #     self.startColumn = 2
    #
    #     # LOX geometry
    #     self.startRow = 1
    #     self.currentColumn = self.startColumn
    #     self.num_elements = 0
    #     while ((self.currentColumn < self.injector_file.ncols) and (
    #             self.injector_file.cell_value(self.startRow, self.currentColumn))):
    #         self.num_elements += 1
    #         self.currentColumn += 1
    #     if (self.num_elements > 0):
    #         self.LOX = LOX()
    #         self.LOX.num_elements = self.num_elements
    #         self.LOX.elements = self.parse_elements(self.LOX.num_elements, self.startRow, self.startColumn)
    #
    #     # Primary Fuel geometry
    #     self.startRow = 8
    #     self.currentColumn = self.startColumn
    #     self.num_elements = 0
    #     while ((self.currentColumn < self.injector_file.ncols) and (
    #             self.injector_file.cell_value(self.startRow, self.currentColumn))):
    #         self.num_elements += 1
    #         self.currentColumn += 1
    #     if (self.num_elements > 0):
    #         self.Fuel = Fuel()
    #         self.Fuel.num_elements = self.num_elements
    #         self.Fuel.elements = self.parse_elements(self.Fuel.num_elements, self.startRow, self.startColumn)
    #
    #     # Secondary Fuel geometry
    #     self.startRow = 15
    #     self.currentColumn = self.startColumn
    #     self.num_elements = 0
    #     while ((self.currentColumn < self.injector_file.ncols) and (
    #             self.injector_file.cell_value(self.startRow, self.currentColumn))):
    #         self.num_elements += 1
    #         self.currentColumn += 1
    #     if (self.num_elements > 0):
    #         self.SecFuel = SecFuel()
    #         self.SecFuel.num_elements = self.num_elements
    #         self.SecFuel.elements = self.parse_elements(self.SecFuel.num_elements, self.startRow, self.startColumn)
    #
    #     #  Window Cooling geometry
    #     self.startRow = 22
    #     self.currentColumn = self.startColumn
    #     self.num_elements = 0
    #     while ((self.currentColumn < self.injector_file.ncols) and (
    #             self.injector_file.cell_value(self.startRow, self.currentColumn))):
    #         self.num_elements += 1
    #         self.currentColumn += 1
    #     if (self.num_elements > 0):
    #         self.WinCool = WinCool()
    #         self.WinCool.num_elements = self.num_elements
    #         self.WinCool.elements = self.parse_elements(self.WinCool.num_elements, self.startRow, self.startColumn)

    # def parse_elements(self, numOfElements, startRow, startColumn):
    #     elements = []
    #     for i in range(numOfElements):
    #         if (self.injector_file.cell_value(startRow, startColumn + i) == "Manifold"):
    #             newElement = Element()
    #             newElement.type = "Manifold"
    #             newElement.length = float('nan')
    #             newElement.diameter = float('nan')
    #             newElement.area = float('nan')
    #             newElement.volume = self.injector_file.cell_value(startRow + 4, startColumn + i)
    #             newElement.numOfElements = self.injector_file.cell_value(startRow + 5, startColumn + i)
    #             newElement.dp = 0
    #         elif (self.injector_file.cell_value(startRow, startColumn + i) == "Orifice"):
    #             newElement = Element()
    #             newElement.type = "Orifice"
    #             newElement.length = float('nan')
    #             newElement.diameter = float('nan')
    #             newElement.area = float('nan')
    #             newElement.volume = float('nan')
    #             newElement.numOfElements = self.injector_file.cell_value(startRow + 5, startColumn + i)
    #             newElement.dp = 0
    #         elif (self.injector_file.cell_value(startRow, startColumn + i) == "Line"):
    #             newElement = Element()
    #             newElement.type = "Line"
    #             newElement.length = self.injector_file.cell_value(startRow + 1, startColumn + i)
    #             newElement.diameter = self.injector_file.cell_value(startRow + 2, startColumn + i)
    #             newElement.area = self.injector_file.cell_value(startRow + 3, startColumn + i)
    #             newElement.volume = float('nan')
    #             newElement.numOfElements = self.injector_file.cell_value(startRow + 5, startColumn + i)
    #             newElement.dp = 0
    #         else:
    #             print("ERROR: Wrong input in Injector Geometry")
    #             sys.exit()
    #
    #         elements.append(newElement)
    #     return elements

    def reduceInjectorHeadToSingleInjector(self):

        # LOX side
        for elementNo in range(len(self.LOX.elements)):
            if (self.LOX.elements[elementNo].type == "Manifold"):
                self.LOX.elements[elementNo].volume = self.LOX.elements[elementNo].volume * self.LOX.elements[
                    elementNo + 1].numOfElements
                self.LOX.elements[elementNo].numOfElements = 1
            elif (self.LOX.elements[elementNo].type == "Orifice"):
                self.LOX.elements[elementNo].numOfElements = 1
            elif (self.LOX.elements[elementNo].type == "Line"):
                self.LOX.elements[elementNo].diameter *= np.sqrt(self.LOX.elements[elementNo].numOfElements)
                self.LOX.elements[elementNo].area = self.LOX.elements[elementNo].diameter ** 2 * np.pi / 4
                self.LOX.elements[elementNo].numOfElements = 1

        # Fuel side
        for elementNo in range(len(self.Fuel.elements)):
            if (self.Fuel.elements[elementNo].type == "Manifold"):
                self.Fuel.elements[elementNo].volume = self.Fuel.elements[elementNo].volume * self.Fuel.elements[
                    elementNo + 1].numOfElements
                self.Fuel.elements[elementNo].numOfElements = 1
            elif (self.Fuel.elements[elementNo].type == "Orifice"):
                self.Fuel.elements[elementNo].numOfElements = 1
            elif (self.Fuel.elements[elementNo].type == "Line"):
                self.Fuel.elements[elementNo].diameter *= np.sqrt(self.Fuel.elements[elementNo].numOfElements)
                self.Fuel.elements[elementNo].area = self.Fuel.elements[elementNo].diameter ** 2 * np.pi / 4
                self.Fuel.elements[elementNo].numOfElements = 1

    # def assignFluidTypes(self):
    #     if (self.LOX):
    #         self.LOX.type = "LOX"
    #     else:
    #         print("ERROR: The program needs a LOX injector geometry")
    #         sys.exit()
    #
    #     if (self.Fuel):
    #         self.Fuel.type = self.config.fuel
    #     else:
    #         print("ERROR: The program needs a primary fuel injector geometry")
    #         sys.exit()
    #
    #     if (hasattr(self, 'SecFuel')):
    #         self.SecFuel.type = self.config.fuel
    #
    #     if (hasattr(self, 'WinCool')):
    #         self.WinCool.type = self.config.fluid_wc

    # def assignTemperaturesAndMassflowsToElements(self):
    #
    #     # LOX
    #     for i in self.LOX.elements:
    #         i.T = self.config.T_io
    #         i.mDot = self.config.mDot_O2 / i.numOfElements
    #
    #     # Fuel
    #     for i in self.Fuel.elements:
    #         i.T = self.config.T_if
    #         i.mDot = self.config.mDot_F / i.numOfElements
    #
    #     # Secondary Fuel
    #     if (hasattr(self, 'SecFuel')):
    #         for i in self.SecFuel.elements:
    #             i.T = self.config.T_sf
    #             i.mDot = self.config.mDot_SF / i.numOfElements
    #
    #     # Window Cooling
    #     if (hasattr(self, 'WinCool')):
    #         for i in self.WinCool.elements:
    #             i.T = self.config.T_wc
    #             i.mDot = self.config.mDot_WC / i.numOfElements

    def calculatePressureInElements(self):

        # Oxidizer
        if math.isnan(self.config.p_io):
            print("P_IO NAN")
            self.iterateUpstreamPressureDrops(self.LOX, 0, self.LOX.num_elements - 1, self.config.p_cc)
        else:
            numElemFirstDome = self.findFirstDome(self.LOX)

            self.calculatePressureDropInFeedSystemSide(self.LOX, numElemFirstDome, self.LOX.num_elements - 1,
                                                       self.config.p_io)
            dp_actual = self.config.p_io - self.config.p_cc
            self.LOX.deltaP = dp_actual
            dp_calculated = self.config.p_io - self.LOX.elements[-1].p + self.LOX.elements[-1].dp
            dp_corr_factor = dp_actual / dp_calculated

            self.correctCalculatedPressureDropInFeedSystemSide(self.LOX, numElemFirstDome, self.LOX.num_elements,
                                                               dp_corr_factor)
            self.updateParamtersForCorrectedPressureInFeedSystem(self.LOX, numElemFirstDome, self.LOX.num_elements - 1)
            if numElemFirstDome > 0:
                self.iterateUpstreamPressureDrops(self.LOX, 0, numElemFirstDome, self.config.p_io)
                self.correctCalculatedPressureDropInFeedSystemSideReverse(self.LOX, 0, numElemFirstDome, dp_corr_factor)

        # Fuel
        if math.isnan(self.config.p_if):
            self.iterateUpstreamPressureDrops(self.Fuel, 0, self.Fuel.num_elements - 1, self.config.p_cc)
        else:
            numElemFirstDome = self.findFirstDome(self.Fuel)
            self.calculatePressureDropInFeedSystemSide(self.Fuel, numElemFirstDome, self.Fuel.num_elements - 1,
                                                       self.config.p_if)
            dp_actual = self.config.p_if - self.config.p_cc
            self.Fuel.deltaP = dp_actual
            dp_calculated = self.config.p_if - self.Fuel.elements[-1].p + self.Fuel.elements[-1].dp
            dp_corr_factor = dp_actual / dp_calculated
            self.correctCalculatedPressureDropInFeedSystemSide(self.Fuel, numElemFirstDome, self.Fuel.num_elements,
                                                               dp_corr_factor)
            # for i in self.Fuel.elements:
            #     print(i.dp)
            self.updateParamtersForCorrectedPressureInFeedSystem(self.Fuel, numElemFirstDome,
                                                                 self.Fuel.num_elements - 1)
            if (numElemFirstDome > 0):
                self.iterateUpstreamPressureDrops(self.Fuel, 0, numElemFirstDome, self.config.p_if)
                self.correctCalculatedPressureDropInFeedSystemSideReverse(self.Fuel, 0, numElemFirstDome,
                                                                          dp_corr_factor)

        if (hasattr(self, 'SecFuel')):
            if (math.isnan(self.config.p_sf)):
                self.iterateUpstreamPressureDrops(self.SecFuel, 0, self.SecFuel.num_elements - 1, self.config.p_cc)
            else:
                numElemFirstDome = self.findFirstDome(self.SecFuel)
                self.calculatePressureDropInFeedSystemSide(self.SecFuel, numElemFirstDome,
                                                           self.SecFuel.num_elements - 1, self.config.p_sf)
                dp_actual = self.config.p_sf - self.config.p_cc
                self.SecFuel.deltaP = dp_actual
                dp_calculated = self.config.p_sf - self.SecFuel.elements[-1].p + self.SecFuel.elements[-1].dp
                dp_corr_factor = dp_actual / dp_calculated
                self.correctCalculatedPressureDropInFeedSystemSide(self.SecFuel, numElemFirstDome,
                                                                   self.SecFuel.num_elements, dp_corr_factor)
                self.updateParamtersForCorrectedPressureInFeedSystem(self.SecFuel, numElemFirstDome,
                                                                     self.SecFuel.num_elements - 1)
                if (numElemFirstDome > 0):
                    self.iterateUpstreamPressureDrops(self.SecFuel, 0, numElemFirstDome, self.config.p_sf)
                    self.correctCalculatedPressureDropInFeedSystemSideReverse(self.SecFuel, 0, numElemFirstDome,
                                                                              dp_corr_factor)

        if (hasattr(self, 'WinCool')):
            if (math.isnan(self.config.p_wc)):
                self.iterateUpstreamPressureDrops(self.WinCool, 0, self.WinCool.num_elements - 1, self.config.p_cc)
            else:
                numElemFirstDome = self.findFirstDome(self.WinCool)
                self.calculatePressureDropInFeedSystemSide(self.WinCool, numElemFirstDome,
                                                           self.WinCool.num_elements - 1, self.config.p_wc)
                dp_actual = self.config.p_wc - self.config.p_cc
                self.WinCool.deltaP = dp_actual
                dp_calculated = self.config.p_wc - self.WinCool.elements[-1].p + self.WinCool.elements[-1].dp
                dp_corr_factor = dp_actual / dp_calculated
                self.correctCalculatedPressureDropInFeedSystemSide(self.WinCool, numElemFirstDome,
                                                                   self.WinCool.num_elements, dp_corr_factor)
                self.updateParamtersForCorrectedPressureInFeedSystem(self.WinCool, numElemFirstDome,
                                                                     self.WinCool.num_elements - 1)
                if (numElemFirstDome > 0):
                    self.iterateUpstreamPressureDrops(self.WinCool, 0, numElemFirstDome, self.config.p_wc)
                    self.correctCalculatedPressureDropInFeedSystemSideReverse(self.WinCool, 0, numElemFirstDome,
                                                                              dp_corr_factor)

    def iterateUpstreamPressureDrops(self, fluid, firstElement, lastElement, p_lastElement):

        upperBoundary = p_lastElement * 3
        lowerBoundary = p_lastElement
        delta = 1e12
        counter = 0
        while (abs(delta) > 0.1):
            if (counter > 200):
                print("Not converged")
                sys.exit()
            counter += 1
            fluid.elements[firstElement].p = (upperBoundary + lowerBoundary) / 2
            self.calculatePressureDropInFeedSystemSide(fluid, firstElement, lastElement, fluid.elements[firstElement].p)
            delta = (fluid.elements[lastElement].p - fluid.elements[lastElement].dp) - p_lastElement
            if (delta > 0):
                upperBoundary = fluid.elements[firstElement].p
            else:
                lowerBoundary = fluid.elements[firstElement].p

    def calculatePressureDropInFeedSystemSide(self, fluid, firstElement, lastElement, p_firstElement):

        if (fluid.type == "CH4"):
            fld = "Methane"
        elif (fluid.type == "H2"):
            fld = "Hydrogen"
        else:
            fld = "Oxygen"

        fluid.elements[firstElement].p = p_firstElement
        fluid.elements[firstElement].rho = getRefPropValue('D', 'T', fluid.elements[firstElement].T, 'P',
                                                           fluid.elements[firstElement].p * 100000, fld)

        if fluid.elements[firstElement].type == "Line":

            fluid.elements[firstElement].velocity = fluid.elements[firstElement].mDot / (
                    fluid.elements[firstElement].rho * fluid.elements[firstElement].area)

            fluid.elements[firstElement].nu = getRefPropValue('V', 'T', fluid.elements[firstElement].T, 'P',
                                                              fluid.elements[firstElement].p * 100000, fld) / \
                                              fluid.elements[firstElement].rho

            fluid.elements[firstElement].dp = func.plinedrop(fluid.elements[firstElement].length,
                                                             fluid.elements[firstElement].diameter,
                                                             fluid.elements[firstElement].rho.fluid.elements[
                                                                 firstElement].velocity,
                                                             fluid.elements[firstElement].nu) * 1e-5

        elif fluid.elements[firstElement].type == "Manifold":
            fluid.elements[firstElement].dp = 0

        for i in range(firstElement + 1, lastElement + 1):
            fluid.elements[i].p = fluid.elements[i - 1].p - fluid.elements[i - 1].dp
            fluid.elements[i].rho = getRefPropValue('D', 'T', fluid.elements[i].T, 'P', fluid.elements[i].p * 100000,
                                                    fld)

            if (i == lastElement):
                if (fluid.elements[i].type == "Orifice"):
                    zeta = 1  # Assuming an instananeous area change with f2/f1 >> 1 and 180°
                    # Pressure for density calculation in last orifice is assumed as CC pressure
                    density = getRefPropValue('D', 'T', fluid.elements[i].T, 'P', self.config.p_cc * 100000, fld)
                    fluid.elements[i].velocity = fluid.elements[i].mDot / (density * fluid.elements[i - 1].area)
                    # fluid.elements[i].velocity = fluid.elements[i].mDot / (fluid.elements[i].rho * fluid.elements[i-1].area)

                    # Old code to calculate the pressure drop
                    fluid.elements[i].dp = func.pdrop(fluid.elements[i].rho, fluid.elements[i].velocity, zeta) * 1e-5
                    # End old approach

                    # New approach to calculate the pressure drop
                    # Respective to the old code, use the chamber pressure p_cc to calculate the density rho
                    # rho = getRefPropValue('D', 'T', fluid.elements[i].T, 'P', self.config.p_cc * 100000, fld)
                    # # rho = fluid.elements[i].rho
                    # C_d = 0.75
                    # # A = fluid.elements[i + 1].area  # line element after orifice
                    # A = fluid.elements[i].area #
                    # mdot = fluid.elements[i].mDot  # mass flow through the orifice
                    # fluid.elements[i].dp = pow(mdot / (C_d * A), 2) / (2 * rho) * 1e-5  # [bar]
                    # End new approach

                else:
                    # Why? -> Recess
                    print("Last element in ", fluid.type, " feed system should be an orifice")

            else:
                if (fluid.elements[i].type == "Orifice"):
                    if ((fluid.elements[i - 1].type == "Manifold") and fluid.elements[i + 1].type == "Line"):
                        zeta = 0.5  # Assuming an instananeous area change with f2/f1 << 1
                        delta = 1e12
                        fluid.elements[i + 1].velocity = fluid.elements[i + 1].mDot / (
                                fluid.elements[i].rho * fluid.elements[i + 1].area)
                        initVelocVal = fluid.elements[i + 1].velocity
                        initRhoVal = fluid.elements[i].rho

                        # New idea to calculate the pressure drop
                        # C_d = 0.75
                        # A = fluid.elements[i+1].area  # line element after orifice
                        # mdot = fluid.elements[i+1].mDot  # mass flow in the line after the orifice equals mass flow through orifice
                        # rho = fluid.elements[i].rho
                        #
                        # fluid.elements[i].dp = pow(mdot / (C_d * A), 2) / (2 * rho) * 1e-5  # [bar]
                        # fluid.elements[i+1].p = fluid.elements[i].p - fluid.elements[i].dp
                        # fluid.elements[i + 1].rho = getRefPropValue('D', 'T', fluid.elements[i + 1].T, 'P',
                        #                                             fluid.elements[i + 1].p * 100000, fld)
                        # End new idea

                        # Old code to calculate the pressure drop
                        try:
                            while (delta > 0.03):
                                old_velocity = fluid.elements[i + 1].velocity
                                fluid.elements[i].dp = func.pdrop(fluid.elements[i].rho, old_velocity, zeta) * 1e-5
                                fluid.elements[i + 1].p = fluid.elements[i].p - fluid.elements[i].dp
                                fluid.elements[i + 1].rho = getRefPropValue('D', 'T', fluid.elements[i + 1].T, 'P',
                                                                            fluid.elements[i + 1].p * 100000, fld)
                                fluid.elements[i].rho = fluid.elements[i + 1].rho
                                fluid.elements[i + 1].velocity = fluid.elements[i + 1].mDot / (
                                        fluid.elements[i + 1].rho * fluid.elements[i + 1].area)
                                delta = abs((fluid.elements[i + 1].velocity - old_velocity) / old_velocity)  # use abs
                        except:
                            fluid.elements[i].dp = func.pdrop(initRhoVal, initVelocVal, zeta) * 1e-5
                        #
                        # End old code

                        fluid.elements[i].velocity = fluid.elements[i + 1].velocity

                    elif ((fluid.elements[i - 1].type == "Line") and (fluid.elements[i + 1].type == "Line")):
                        angle = math.pi / 2  # Assuming 180° angle at cross section
                        zeta = func.zeta_p(angle, fluid.elements[i - 1].area, fluid.elements[i + 1].area)
                        fluid.elements[i].velocity = fluid.elements[i].mDot / (
                                    fluid.elements[i].rho * fluid.elements[i - 1].area)  # v = m / (rho * A_{-1})

                        # Old code to calculate the pressure drop
                        fluid.elements[i].dp = func.pdrop(fluid.elements[i].rho, fluid.elements[i].velocity,
                                                          zeta) * 1e-5
                        # End old code

                        #
                        # New idea to calculate the pressure drop
                        # C_d = 0.75
                        # A = fluid.elements[i + 1].area  # line element after orifice
                        # mdot = fluid.elements[i + 1].mDot
                        # rho = fluid.elements[i].rho
                        # fluid.elements[i].dp = pow(mdot / (C_d * A), 2) / (2 * rho) * 1e-5
                        #
                        # End new idea

                    elif ((fluid.elements[i - 1].type == "Line") and (fluid.elements[i + 1].type == "Manifold")):
                        zeta = 1  # Assuming instantaneous area change with f2/f1 >> 1 and 180°
                        fluid.elements[i].velocity = fluid.elements[i].mDot / (
                                fluid.elements[i].rho * fluid.elements[i - 1].area)
                        fluid.elements[i].dp = func.pdrop(fluid.elements[i].rho, fluid.elements[i].velocity,
                                                          zeta) * 1e-5

                    else:
                        print("Orifice in element ", i, " in ", fluid.type,
                              " feed system must be surrounded by Line/Mainfold")

                elif (fluid.elements[i].type == "Line"):
                    fluid.elements[i].velocity = fluid.elements[i].mDot / (
                                fluid.elements[i].rho * fluid.elements[i].area)
                    fluid.elements[i].nu = getRefPropValue('V', 'T', fluid.elements[i].T, 'P',
                                                           fluid.elements[i].p * 100000, fld) / fluid.elements[i].rho

                    fluid.elements[i].dp = func.plinedrop(fluid.elements[i].length, fluid.elements[i].diameter,
                                                          fluid.elements[i].rho, fluid.elements[i].velocity,
                                                          fluid.elements[i].nu) * 1e-5

                elif (fluid.elements[i].type == "Manifold"):
                    fluid.elements[i].dp = 0

    def findFirstDome(self, fluid):

        for i in range(fluid.num_elements):
            if (fluid.elements[i].type == "Manifold"):
                numElem = i
                break
            else:
                numElem = float('nan')

        return numElem

    def correctCalculatedPressureDropInFeedSystemSide(self, fluid, firstElement, lastElement, dp_corr_factor):

        fluid.elements[firstElement].dp *= dp_corr_factor

        for i in range(firstElement + 1, lastElement):
            fluid.elements[i].p = fluid.elements[i - 1].p - fluid.elements[i - 1].dp
            fluid.elements[i].dp *= dp_corr_factor

    def updateParamtersForCorrectedPressureInFeedSystem(self, fluid, firstElement, lastElement):

        if (fluid.type == "CH4"):
            fld = "Methane"
        elif (fluid.type == "H2"):
            fld = "Hydrogen"
        else:
            fld = "Oxygen"

        fluid.elements[firstElement].rho = getRefPropValue('D', 'T', fluid.elements[firstElement].T, 'P',
                                                           fluid.elements[firstElement].p * 100000, fld)

        if (fluid.elements[firstElement].type == "Line"):
            fluid.elements[firstElement].velocity = fluid.elements[firstElement].mDot / (
                    fluid.elements[firstElement].rho * fluid.elements[firstElement].area)
            fluid.elements[firstElement].nu = getRefPropValue('V', 'T', fluid.elements[firstElement].T, 'P',
                                                              fluid.elements[firstElement].p * 100000, fld) / \
                                              fluid.elements[firstElement].rho

        for i in range(firstElement + 1, lastElement + 1):
            fluid.elements[i].rho = getRefPropValue('D', 'T', fluid.elements[i].T, 'P', fluid.elements[i].p * 100000,
                                                    fld)

            if (i == lastElement):
                if (fluid.elements[i].type == "Orifice"):
                    fluid.elements[i].velocity = fluid.elements[i].mDot / (
                            fluid.elements[i].rho * fluid.elements[i - 1].area)
                else:
                    print("Last element in ", fluid.type, " feed system should be an orifice")
            else:
                if (fluid.elements[i].type == "Orifice"):
                    if ((fluid.elements[i - 1].type == "Manifold") and fluid.elements[i + 1].type == "Line"):
                        fluid.elements[i + 1].velocity = fluid.elements[i + 1].mDot / (
                                fluid.elements[i].rho * fluid.elements[i + 1].area)
                        fluid.elements[i].velocity = fluid.elements[i + 1].velocity
                    elif ((fluid.elements[i - 1].type == "Line") and (fluid.elements[i + 1].type == "Line")):
                        fluid.elements[i].velocity = fluid.elements[i].mDot / (
                                fluid.elements[i].rho * fluid.elements[i - 1].area)
                    elif ((fluid.elements[i - 1].type == "Line") and (fluid.elements[i + 1].type == "Manifold")):
                        fluid.elements[i].velocity = fluid.elements[i].mDot / (
                                fluid.elements[i].rho * fluid.elements[i - 1].area)
                    else:
                        print("Orifice in element ", i, " in ", fluid.type,
                              " feed system must be surrounded by Line/Mainfold")
                elif (fluid.elements[i].type == "Line"):
                    fluid.elements[i].velocity = fluid.elements[i].mDot / (
                            fluid.elements[i].rho * fluid.elements[i].area)
                    fluid.elements[i].nu = getRefPropValue('V', 'T', fluid.elements[i].T, 'P',
                                                           fluid.elements[i].p * 100000, fld) / fluid.elements[i].rho

        # print("PressureDrop-Dist")
        # for i in fluid.elements:
        #     try:
        #         print(i.dp, i.rho)
        #     except:
        #         print(0)
        # print("")

    def correctCalculatedPressureDropInFeedSystemSideReverse(self, fluid, firstElement, lastElement, dp_corr_factor):

        for i in range(lastElement - 1, firstElement, -1):
            fluid.elements[i].dp *= dp_corr_factor
            fluid.elements[i].p = fluid.elements[i + 1].p + fluid.elements[i].dp

    def calculateAcousticParameters(self):

        # LOX
        for i in range(self.LOX.num_elements):
            self.updateAcousticParametersElement(self.LOX.elements[i], self.LOX.type, self.config.mDot_O2,
                                                 self.config.p_cc)

        # Fuel
        for i in range(self.Fuel.num_elements):
            self.updateAcousticParametersElement(self.Fuel.elements[i], self.Fuel.type, self.config.mDot_F,
                                                 self.config.p_cc)

        # SecFuel
        if (hasattr(self, 'SecFuel')):
            for i in range(self.SecFuel.num_elements):
                self.updateAcousticParametersElement(self.SecFuel.elements[i], self.SecFuel.type, self.config.mDot_SF,
                                                     self.config.p_cc)

        # WinCool
        if (hasattr(self, 'WinCool')):
            for i in range(self.WinCool.num_elements):
                self.updateAcousticParametersElement(self.WinCool.elements[i], self.WinCool.type, self.config.mDot_WC,
                                                     self.config.p_cc)

    def updateAcousticParametersElement(self, element, fluid_type, mDot_SystemSide, p_cc):
        ## Not used when the flow solver is turned on !!!!!!!

        if (fluid_type == "CH4"):
            fld = "Methane"
        elif (fluid_type == "H2"):
            fld = "Hydrogen"
        else:
            fld = "Oxygen"

        if (element.type == "Line"):
            element.rho = getRefPropValue('D', 'T', element.T, 'P', element.p * 100000, fld)
            #element.velocity = element.mDot / (element.rho * element.area)
            #element.soundVelocity = getRefPropValue('A', 'T', element.T, 'P', element.p * 100000, fld)
            #element.t = element.length / element.soundVelocity
            #element.Z = element.soundVelocity * mDot_SystemSide / (element.area * p_cc * 1e5)
            #element.M = element.velocity / element.soundVelocity

        elif (element.type == "Orifice"):
            element.rho = getRefPropValue('D', 'T', element.T, 'P', element.p * 100000, fld)
            #element.soundVelocity = getRefPropValue('A', 'T', element.T, 'P', element.p * 100000, fld)
            #element.Z = (2 * element.dp * 1e5 / element.mDot) * (mDot_SystemSide / p_cc * 1e5)

        elif (element.type == "Manifold"):
            element.rho = getRefPropValue('D', 'T', element.T, 'P', element.p * 100000, fld)
            #element.soundVelocity = getRefPropValue('A', 'T', element.T, 'P', element.p * 100000, fld)
            #element.K = math.pow(element.soundVelocity, 2) * element.rho
            #element.C = element.rho * element.volume * p_cc * 1e5 / (element.K * element.mDot)

    def getLengthCas(self, fluid):

        accumLength = 0

        if (fluid == "Oxygen"):
            for i in self.LOX.elements:
                if (i.type == "Line"):
                    accumLength += i.length
        elif (fluid == "Fuel"):
            for i in self.Fuel.elements:
                if (i.type == "Line"):
                    accumLength += i.length
        else:
            print("ERROR - Unknown fluid")
            sys.exit()

        return accumLength

#     def calculate_feed_systems_parameter(self):
#
#         from sympy import sqrt, log
#         from scipy.optimize import newton_krylov
#
#         # New approach to model the feed system, this includes
#         # - the propellant tank
#         # - the connection(s) between the propellant tank and the injector
#         # - the injector itself
#         #
#         # Assumptions:
#         # - The first element is considered to be the propellant/oxidizer tank (and therefore be a manifold)
#         # - The temperature is constant between the tank and the start of the injector (T_t), and further constant over
#         #   the whole injector (T_i)
#         # - No secondary fuel line and no window cooling option are considered (to be done, if the results improve)
#         # - The oxidizer/fuel combination is LOX/Hydrogen
#         #
#         # Warning - the calculation using the whole feed system will need more time (Don't know why yet).
#         #
#         # To use this, it is necessary to
#         # - provide values for p_ot, p_ft, T_ot, P_ft in the settings file
#         # - set the use_whole_feed_system parameter in the setting file to True. If this parameter is set to False, the
#         #   tool will use the 'previous' code to calculate the injector parameters (see line 150)
#         # - modify the injector geometry file to include additional elements for the feed system
#         # - specify the first element of the actual injector in the settings file
#         #
#         # Defining the symbols used here
#         #   d           diameter, or hydraulic diameter. For a cylindrical cross section, the hydraulic diameter and the
#         #               actual diameter are equal. For a concentric ring gap, the hydraulic diameter is d = d_out - d_in
#         #   l           length of a pipe (line) element [m]
#         #   A           cross section area of an element [m^2]
#         #   p           local static pressure [bar]
#         #   dp          local static pressure drop over an element [bar]
#         #   p_total     total pressure [bar]
#         #   dp_total    local total pressure loss over an element [bar]
#         dp_total = lambda zeta, rho, v: zeta * rho * pow(v, 2) / 2
#         #   v           local velocity [m/s]
#         #   rho         local density [kg/m^3]
#         #   eta         dynamic viscosity [kg/(m s)]
#         #   Re          local reynolds number [-]
#         Re = lambda rho, v, d, eta: rho * v * d / eta
#         #
#         #   epsilon     coefficient of friction [-], dependent on geometry and flow
#         #               Assuming linear flow for Re < 2300, epsilon is only dependent on Re,
#         # epsilon_zylinder = lambda Re: 64 / Re
#         # epsilon_concentric_ring_gap = lambda Re: 96 / Re
#         #                For a turbulent flow Re > 2300, epsilon can be calculated using the Colebrook formula. This
#         #                has to be solved numerically, here, the newton_krylov approach with an initial guess of is
#         #                used. The initial guess should not affect the results. (log(x, 10) is the logarithm to base 10)
#         epsilon = lambda d, Re, kw: (float(newton_krylov(
#             lambda epsilon: 1 / sqrt(epsilon) + 2 * log(2.51 / (sqrt(epsilon) * Re), 10) - 0.27 * d / kw, 0.25)))
#         #
#         #   zeta        pressure loss factor for the total pressure  [-]
#         zeta = lambda epsilon, l, d: epsilon * l / d
#         #   c           local speed of sound [m/s]
#         #   Ma          Mach number [-]
#         Ma = lambda v, c: v / c
#         #
#         #   k           adiabatic index [-]
#         #   T_total     total temperature [K]
#         #   T_adiabatic local temperature assuming an adiabatic process [K]
#         T_adiabatic = lambda T_total, k, Ma: T_total / (1 + (k - 1) / 2 * pow(Ma, 2))
#         #   p_isentropic    local static pressure calculated assuming an isentropic process
#         p_isentropic = lambda p_total, k, Ma: p_total / ((1 + (k - 1) / 2 * pow(Ma, 2)) ^ (k / (k - 1)))
#
#         # Performing the calculation for the LOX part, and then for the Fuel part of the feed system.
#         # For more simplicity, a secondary fuel part, or a window part are not considered in the first approach here.
#         for fluid_part in [self.LOX, self.Fuel]:
#
#             # For more simplicity, define the known values as new variables
#             elements = fluid_part.elements  # Just for better readability (the name of the variable is now shorter)
#             p_cc = self.config.p_cc  # Pressure in the combustion chamber [bar]
#             if fluid_part == self.LOX:
#                 fluid = 'Oxygen'
#                 p_t = self.config.p_ot  # Pressure in the oxidizer tank [bar]
#                 T_t = self.config.T_ot  # Temperature in the oxidizer tank [K]
#                 p_i = self.config.p_io  # Pressure in the first manifold of the actual injector [bar]
#                 T_i = self.config.T_io  # Temperature in the first manifold of the actual injector [K]
#                 mdot = self.config.mDot_O2  # Mass flow for the whole oxidizer feed system [kg/s]
#             elif fluid_part == self.Fuel:
#                 fluid = 'Hydrogen'
#                 p_t = self.config.p_ft  # Pressure in the oxidizer tank [bar]
#                 T_t = self.config.T_ft  # Temperature in the oxidizer tank [K]
#                 p_i = self.config.p_if  # Pressure in the first manifold of the actual injector [bar]
#                 T_i = self.config.T_if  # Temperature in the first manifold of the actual injector [K]
#                 mdot = self.config.mDot_F  # Mass flow for the whole fuel feed system [kg/s]
#
#             # Finding the injector start
#             # This information will be used to split the feed system into everything before the injector,
#             # and the actual injector.
#             injector_start_index = self.config.injector_start_index
#
#             # Check if there are any elements before the actual injector. If not, use the injector values for the
#             # tank values
#             if injector_start_index == 0:
#                 p_t = p_i
#                 T_t = T_i
#
#             # Propellant/Oxidizer tank The first element (equals index zero = [0]) in the given fluid_part is
#             # considered to be the propellant/oxidizer tank, and is therefore a manifold with known parameters.
#             elements[0].p = p_t  # Absolute pressure of the propellant tank [bar]
#             elements[0].dp = 0  # No pressure drop in/over the propellant/oxidizer tank is considered [bar]
#             elements[0].T = T_t  # Temperature in the propellant/oxidizer tank [K]
#             elements[0].rho = getRefPropValue('D', 'T', T_t, 'P', p_t * 100000, fluid)  # Density [kg/m^3]
#
#             # As the element is a manifold, there's no velocity. Therefore, the total values euqal the static values.
#             elements[0].velocity = 0  # [m/s]
#             elements[0].p_total = p_t  # [bar]
#             elements[0].T_total = T_t  # [K]
#             elements[0].rho_total = elements[0].rho  # [kg/m^3]
#             elements[0].soundVelocity = getRefPropValue('A', 'T', element.T, 'P', element.p * 100000, fluid)  # [m/s]
#
#             # Iterating over all elements of the current fluid part, without including the first element - which is
#             # the propellant/oxidizer manifold.
#             # This calculates the pressure drop over each element (and all necessary parameters to do so)
#             for i in range(1, len(elements)):
#
#                 # Absolute pressure
#                 if i == injector_start_index:
#                     # If the current element is the start of the injector, the pressure is provided by the settings file
#                     elements[i].p = p_i
#                 else:
#                     # Based on the absolute pressure and the pressure drop of the previous element, the pressure of the
#                     # current element can be calculated
#                     elements[i].p = elements[i - 1].p - elements[i - 1].dp
#
#                 # If the current element is an element before the injector, the temperature of this element has to be
#                 # updated. If the element is within the actual injector, the temperature has already been set
#                 # accordingly (see line 121 - 128)
#                 if i < injector_start_index:
#                     elements[i].T = T_t
#
#                 # As pressure and temperature are known, the density is determined
#                 elements[i].rho = getRefPropValue('D', 'T', elements[i].T, 'P', elements[i].p * 100000, fluid)
#
#                 if elements[i].type == 'Line':
#                     # Velocity is calculated using the constant mass flow m = rho * A * v = const. [m/s]
#                     elements[i].velocity = elements[i].mDot / (elements[i].rho * elements[i].area)
#                     # Kinematic viscosity by dividing the dynamic viscosity by the density
#                     elements[i].nu = getRefPropValue('V', 'T', elements[i].T, 'P', elements[i].p * 100000, fluid) / \
#                                      elements[i].rho
#                     # Pressure drop over the line element [bar]
#                     elements[i].dp = func.plinedrop(elements[i].length, elements[i].diameter, elements[i].rho,
#                                                     elements[i].velocity, elements[i].nu) * 1e-5
#
#                 elif elements[i].type == 'Manifold':
#                     elements[i].dp = 0  # No pressure drop on a manifold [bar]
#
#                 elif elements[i].type == 'Orifice':
#                     # Multiple approaches to calculate the pressure drop over an orifice are possible.
#                     # One has to be selected.
#
#                     # Approach No. 1
#                     # This would require to assign a cross section area to an orifice element in the config file
#                     # dp = (mdot / (C_d * A))^2 / (2 * rho)
#                     #
#                     # C_d = 0.75
#                     # elements[i].dp = pow(elements[i].mDot / (C_d * elements[i].area), 2) / (2 * elements[i].rho) * 1e-5  # [bar]
#
#                     # Approach No. 2
#                     # This is the approach used by Geuking in the previous version.
#                     # The calculation/formula is based on the different combination of element, possibilities are
#                     # - Line, Orifice -> End of the injector. If this is the case, then i == len(elements) - 1
#                     # - Line, Orifice, Line
#                     # - Line, Orifice, Manifold
#                     # - Manifold, Orifice, Line
#
#                     if i == len(elements) - 1 or (
#                             elements[i - 1].type == 'Line' and elements[i + 1].type == 'Manifold'):
#                         # - Line, Orifice (end of injector)
#                         # - Line, Orifice, Manifold
#                         #
#                         # Using zeta as a loss factor based on the geometry, assuming an instantaneous area change with
#                         # 180° and A[i-1] / A[i] >> 1
#                         zeta = 1
#                         elements[i].velocity = elements[i].mDot / (elements[i].rho * elements[i - 1].area)
#                         elements[i].dp = func.pdrop(elements[i].rho, elements[i].velocity, zeta) * 1e-5  # [bar]
#
#                     elif elements[i - 1].type == 'Line' and elements[i + 1].type == 'Line':
#                         # - Line, Orifice, Line
#                         #
#                         # Using zeta as a loss factor based on the geometry, assuming a 180° angle between the lines
#                         zeta = func.zeta_p(math.pi / 2, elements[i - 1].area, elements[i + 1].area)
#                         elements[i].velocity = elements[i].mDot / (elements[i].rho * elements[i - 1].area)
#                         elements[i].dp = func.pdrop(elements[i].rho, elements[i].velocity, zeta) * 1e-5  # [bar]
#
#                     elif elements[i - 1].type == 'Manifold' and elements[i + 1].type == 'Line':
#                         # - Manifold, Orifice, Line
#                         #
#                         # This iteration is quite dangerous, there's no guarantee this converges at all. In a worst
#                         # case, this is an infinite loop.
#                         #
#                         # According to the formula used here, the pressure drop is dependent on rho[i+1] - which itself
#                         # depends on the pressure drop. The iteration tries to approximates this.
#                         #
#                         # Using zeta as a loss factor based on the geometry, assuming an instantaneous area change with
#                         # 180° and A[i-1] / A[i] << 1
#                         zeta = 0.5
#                         delta = 1e12
#                         elements[i + 1].velocity = elements[i + 1].mDot / (elements[i].rho * elements[i + 1].area)
#                         init_velocity = elements[i + 1].velocity
#                         init_rho = elements[i].rho
#                         while delta > 0.03:
#                             old_velocity = elements[i + 1].velocity
#                             elements[i].dp = func.pdrop(elements[i].rho, old_velocity, zeta) * 1e-5
#                             elements[i + 1].p = elements[i].p - elements[i].dp
#                             elements[i + 1].rho = getRefPropValue('D', 'T', elements[i + 1].T, 'P',
#                                                                   elements[i + 1].p * 100000, fluid)
#                             elements[i].rho = elements[i + 1].rho
#                             elements[i + 1].velocity = elements[i + 1].mDot / (
#                                         elements[i + 1].rho * elements[i + 1].area)
#                             delta = abs((elements[i + 1].velocity - old_velocity) / old_velocity)
#
#                 else:
#                     # Just for extra error checking, this should never occur. Crash on purpose.
#                     print("Error while calculating the injector. There's a wrong element within the injector.")
#                     sys.exit()
#
#             # Create documentation
#             injector_documentation(elements, os.path.join(self.config.result_path, 'Injector'),
#                                    '{}_calculated.csv'.format(fluid), p_t, p_i, p_cc, T_t, T_i)
#
#             # To account for small errors in the pressure drop estimations, the pressure drop is linear 'scaled' to
#             # match p_t, p_i and p_cc.
#             # Scaling all elements before the actual injector - only if there are elements before the injector
#             if injector_start_index != 0:
#                 dp_ideal = p_t - p_i
#                 dp_calculated = elements[0].p - elements[injector_start_index - 1].p + elements[
#                     injector_start_index - 1].dp
#                 correction = dp_ideal / dp_calculated
#                 elements[0:injector_start_index] = pressure_drop_correction(elements[0:injector_start_index],
#                                                                             correction)
#             # Scaling all elements within the actual injector
#             dp_ideal = p_i - p_cc
#             dp_calculated = elements[injector_start_index].p - elements[-1].p + elements[-1].dp
#             correction = dp_ideal / dp_calculated
#             elements[injector_start_index:] = pressure_drop_correction(elements[injector_start_index:], correction)
#
#             # Iterate over all elements and calculate all further necessary parameters which are needed to create the
#             # transfer functions later
#             for element in elements:
#
#                 element.rho = getRefPropValue('D', 'T', element.T, 'P', element.p * 100000, fluid)
#                 element.soundVelocity = getRefPropValue('A', 'T', element.T, 'P', element.p * 100000, fluid)
#
#                 if element.type == 'Line':
#                     element.velocity = element.mDot / (element.rho * element.area)
#                     element.t = element.length / element.soundVelocity
#                     element.Z = element.soundVelocity * mdot / (element.area * p_cc * 1e5)
#                     element.M = element.velocity / element.soundVelocity
#
#                 elif element.type == 'Manifold':
#                     element.K = pow(element.soundVelocity, 2) * element.rho
#                     element.C = element.rho * element.volume * p_cc * 1e5 / (element.K * element.mDot)
#
#                 elif element.type == 'Orifice':
#                     element.Z = (2 * element.dp / element.mDot) * (mdot / p_cc)
#
#             # Create documentation
#             injector_documentation(elements, os.path.join(self.config.result_path, 'Injector'),
#                                    '{}_scaled.csv'.format(fluid), p_t, p_i, p_cc, T_t, T_i)
#
#
# def acoustic_parameters(elements, mdot, p_cc):
#     for element in elements:
#         if element.type == 'Line':
#             element.velocity = element.mDot / (element.rho * element.area)
#             element.t = element.length / element.soundVelocity
#             element.Z = element.soundVelocity * mdot / (element.area * p_cc * 1e5)
#             element.M = element.velocity / element.soundVelocity
#
#         elif element.type == 'Manifold':
#             element.K = pow(element.soundVelocity, 2) * element.rho
#             element.C = element.rho * element.volume * p_cc * 1e5 / (element.K * element.mDot)
#
#         elif element.type == 'Orifice':
#             element.Z = (2 * element.dp / element.mDot) * (mdot / p_cc)
#
#     return elements
#
#
# def pressure_drop_correction(elements, factor):
#     """ Scales the pressure drops, and corrects depending parameters
#
#     :param elements: List of elements for which the pressure drop has to be corrected. After the pressure drop is
#     scaled, other parameters are re-calculated.
#     :param factor: Correction factor for the pressure drop.
#     :return: The same list of elements, but with corrected parameters .
#
#     :type elements: list[Element]
#     :type factor: float
#     :rtype: list[Element]
#     """
#
#     elements[0].dp *= factor
#
#     for i in range(1, len(elements)):
#         elements[i].dp *= factor
#         elements[i].p = elements[i - 1].p - elements[i - 1].dp
#         # elements[i].rho = getRefPropValue('D', 'T', elements[i].T, 'P', elements[i].p * 100000, fluid)
#
#         # Pretty sure these parameters are not necessary even if they were calculated in the previous code
#         #
#         # if elements[i].type == 'Line':
#         #     elements[i].velocity = elements[i].mDot / (elements[i].rho * elements[i].area)
#         #     # Not sure if the parameter nu is needed
#         #     elements[i].nu = getRefPropValue('V', 'T', elements[i].T, 'P', elements[i].p * 100000, fluid)
#         # elif elements[i].type == 'Orifice':
#         #     # Not sure if these parameters are needed. However, the previous code re-calculated them
#         #     if elements[i-1].type == 'Manifold':
#         #         elements[i].velocity = elements[i].mDot / (elements[i].rho * elements[i+1].area)
#         #     else:
#         #         elements[i].velocity = elements[i].mDot / (elements[i].rho * elements[i-1].area)
#
#     return elements
#
#
# def injector_documentation(elements, path, file, p_t, p_i, p_cc, T_t, T_i):
#     """ Creates documentation about the injector.
#
#     :param elements: List of all elements
#     :param path: path to store results
#     :param file: filename
#
#     :type elements: list[Element]
#     :type file: str
#     :rtype: None
#     """
#
#     import pathlib
#     pathlib.Path(path).mkdir(parents=True, exist_ok=True)
#
#     # Collecting all data
#     data = []
#     for element in elements:
#         data.append([
#             getattr(element, 'type', None),
#             getattr(element, 'area', None),
#             getattr(element, 'diameter', None),
#             getattr(element, 'volume', None),
#             getattr(element, 'length', None),
#             getattr(element, 'numOfElements', None),
#             getattr(element, 'rho', None),
#             getattr(element, 'p', None),
#             getattr(element, 'dp', None),
#             getattr(element, 'T', None),
#             getattr(element, 'velocity', None),
#             getattr(element, 'mDot', None)
#         ])
#
#     data = pd.DataFrame(data)
#     data.columns = ['type', 'area', 'diameter', 'volume', 'length', 'number_of_elements', 'rho', 'p', 'dp', 'T',
#                     'velocity', 'mdot']
#
#     # Saving the data for later use
#     data.to_csv(os.path.join(path, file), index=False)
#
#
# def injector_plots(path, file_calculated, file_scaled, file_result, p_t, p_i, p_cc):
#     import matplotlib.pyplot as plt
#
#     # Reading the data
#     data_calculated = pd.read_csv(os.path.join(path, file_calculated))
#     data_scaled = pd.read_csv(os.path.join(path, file_scaled))
#
#     x = [0]
#     for element in data_scaled['type']:
#         if element == 'Manifold':
#             x.append(x[-1] + 1)
#         if element == 'Line':
#             x.append(x[-1] + 1)
#         if element == 'Orifice':
#             x.append(x[-1])
#
#     plt.figure()
#     plt.plot(x, data_calculated['p'].append(pd.Series(data_calculated['p'].iloc[-1] - data_calculated['dp'].iloc[-1]),
#                                             ignore_index=True), label='calculated')
#     plt.plot(x, data_scaled['p'].append(pd.Series(data_scaled['p'].iloc[-1] - data_scaled['dp'].iloc[-1]),
#                                         ignore_index=True), label='scaled')
#     x_min = min(plt.xlim())
#     x_max = max(plt.xlim())
#     plt.hlines(p_t, xmin=x_min, xmax=x_max, colors='black', linestyles=':', alpha=0.5)
#     plt.hlines(p_i, xmin=x_min, xmax=x_max, colors='black', linestyles=':', alpha=0.5)
#     plt.hlines(p_cc, xmin=x_min, xmax=x_max, colors='black', linestyles=':', alpha=0.5)
#     plt.title('static pressure')
#     plt.ylabel('p [bar]')
#     plt.legend()
#     plt.xticks([])
#     plt.savefig(os.path.join(path, file_result))

def new_manifold(number_of_elements, mDot, volume, c, rho, p_end):
    """ Creates a new manifold injector element

    :param number_of_elements: number of elements
    :param mDot: mass flow through the element [kg/s]
    :param volume: manifold volume [m^3]
    :param c: sound velocity [m/s]
    :param rho: density [kg/m^3]
    :param p_end: pressure in the combustion chamber (or after elements) [Pa]
    :return: a manifold element
    :type number_of_elements: int
    :type mDot: float
    :type volume: float
    :type c: float
    :type rho
    :type p_end: float
    :rtype: Element

    """

    # Mass flow through the whole system, independent from the number of elements
    mDot_system = mDot * number_of_elements
    TYPE = 'Manifold'
    data = pd.Series({'type': TYPE, 'num_elem': number_of_elements, 'volume': volume, 'length': None,
                      'hydraulic diameter': None, 'area': None})
    manifold = Element(data, None, mDot_system)
    manifold.rho = rho
    manifold.soundVelocity = c

    # Acoustic Parameters
    manifold.K = pow(manifold.soundVelocity, 2) * manifold.rho
    manifold.C = manifold.rho * manifold.volume * p_end / (manifold.K * mDot)
    manifold.ImpC = manifold.rho * manifold.volume / manifold.K

    return manifold


def new_orifice(number_of_elements, mDot, dp, p_end, v=None):
    """ Create a new orifice injector element

    :param number_of_elements: number of elements
    :param mDot: mass flow through the element [kg/s]
    :param dp: local static pressure drop over the orifice [Pa]
    :param p_end: pressure in the combustion chamber (or after elements) [Pa]
    :param v: velocity at the end of the orifice. Only needs to be provided if the orifice is at the end of the system
    :return: a orifice element
    :type number_of_elements: int
    :type mDot: float
    :type dp: float
    :type p_end: float
    :rtype: Element

    """

    # Mass flow through the whole system, independent from the number of elements
    mDot_system = mDot * number_of_elements
    TYPE = 'Orifice'
    data = pd.Series({'type': TYPE, 'num_elem': number_of_elements, 'volume': None, 'length': None,
                      'hydraulic diameter': None, 'area': None})
    orifice = Element(data, None, mDot_system)
    orifice.dp = dp

    # Acoustic parameters
    # element.Z = (2 * element.dp / element.mDot) * (mdot / p_cc)
    orifice.Z = (2 * orifice.dp / orifice.mDot) * (mDot_system / p_end)
    orifice.R0 = 2 * orifice.dp / orifice.mDot

    if v is not None:
        orifice.velocity = v

    return orifice


def new_line(number_of_elements, mDot, A, v, c, length, p_end):
    """ Creates a new line injector element

    :param number_of_elements: number of elements
    :param mDot: mass flow through the element [kg/s]
    :param A: cross section area [m^2]
    :param v: velocity [m/s]
    :param c: sound velocity [m/s]
    :param length: length of the line [m]
    :param p_end: pressure in the combustion chamber (or after elements) [Pa]
    :return: a orifice element
    :rtype: Element
    :type number_of_elements: int
    :type mDot: float
    :type A: float
    :type v: float
    :type c: float
    :type length: float
    :type p_end: float
    :rtype: Element
    """

    # Mass flow through the whole system, independent from the number of elements
    mDot_system = mDot * number_of_elements

    TYPE = 'Line'
    data = pd.Series({'type': TYPE, 'num_elem': number_of_elements, 'volume': None, 'length': length,
                      'hydraulic diameter': None, 'area': A})
    line = Element(data, None, mDot_system)
    line.soundVelocity = c
    line.velocity = v

    # Acoustic parameters
    line.t = line.length / line.soundVelocity
    line.Z = line.soundVelocity / line.area * (mDot_system / p_end)
    line.M = line.velocity / line.soundVelocity
    line.dp = None
    line.L0 = line.length / line.area
    return line


def save_injector_elements(path, name, elements):
    """

    """

    with open(os.path.join(path, '{}.txt'.format(name)), 'w') as file:
        file.write(
            '{:<13}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<17}{:<10}{:<10}{:<10}{:<15}{:<15}\n'.format('type', 'noe', 'volume',
                                                                                                'length', 'area',
                                                                                                'massflow', 'dp',
                                                                                                'velocity',
                                                                                                'rho',
                                                                                                'sound_velocity', 't * 1e-4',
                                                                                                'Z', 'M', 'C * 1e-7', 'K'))

        for element in elements:
            element_type = getattr(element, 'type', None)
            noe = getattr(element, 'numOfElements', None)
            volume = getattr(element, 'volume', None)
            if volume:
                volume = np.round(volume * 1e9, 4)
            length = getattr(element, 'length', None)
            if length:
                length = np.round(length * 1e3, 4)
            area = getattr(element, 'area', None)
            if area:
                area = np.round(area * 1e6, 4)
            massflow = getattr(element, 'mDot', None)
            if massflow:
                massflow = np.round(massflow, 4)
            dp = getattr(element, 'dp', None)
            if dp:
                dp = np.round(dp, 4)
            velocity = getattr(element, 'velocity', None)
            if velocity:
                velocity = np.round(velocity, 4)
            rho = getattr(element, 'rho', None)
            if rho:
                rho = np.round(rho, 4)
            sound_velocity = getattr(element, 'soundVelocity', None)
            if sound_velocity:
                sound_velocity = np.round(sound_velocity, 4)
            t = getattr(element, 't', None)
            if t:
                t = np.round(t * 1e4, 4)
            Z = getattr(element, 'Z', None)
            if Z:
                Z = np.round(Z, 4)
            M = getattr(element, 'M', None)
            if M:
                M = np.round(M, 4)
            C = getattr(element, 'C', None)
            if C:
                C = np.round(C * 1e7, 4)
            K = getattr(element, 'K', None)
            if K:
                K = np.round(K, 4)

            file.write(
                '{:<13}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<17}{:<10}{:<10}{:<10}{:<15}{:<15}\n'.format(str(element_type),
                                                                                                    str(noe),
                                                                                                    str(volume),
                                                                                                    str(length),
                                                                                                    str(area),
                                                                                                    str(massflow),
                                                                                                    str(dp),
                                                                                                    str(velocity),
                                                                                                    str(rho),
                                                                                                    str(sound_velocity),
                                                                                                    str(t), str(Z),
                                                                                                    str(M), str(C),
                                                                                                    str(K)))
