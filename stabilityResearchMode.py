# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 15:58:46 2020

@author: micha
"""
import sys, os
import ConcentratedCombustionModel as CCM
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import xlrd
import copy

from matplotlib import cm
from Injector import injector_configuration


# class temperatureVariation:
#
#     def __init__(self, config):
#         self.newConfig = config
#
#         self.run()
#         self.createTempVariationPlots()
#
#     def run(self):
#         # Define temperature range
#         startTempOx = 0.8 * self.newConfig.T_io
#         endTempOx = 1.2 * self.newConfig.T_io
#         self.tempRangeOx = np.linspace(startTempOx, endTempOx, 10)
#
#         startTempFu = 0.8 * self.newConfig.T_if
#         endTempFu = 1.2 * self.newConfig.T_if
#         self.tempRangeFu = np.linspace(startTempFu, endTempFu, 10)
#
#         # Setting parameter values as average of input values
#         self.newConfig.tau_Ox = [(self.newConfig.tau_Ox[0] + self.newConfig.tau_Ox[1]) / 2,
#                                  (self.newConfig.tau_Ox[0] + self.newConfig.tau_Ox[1]) / 2]
#         self.newConfig.dtau_Ox = 1
#         self.newConfig.tau_Fu = [(self.newConfig.tau_Fu[0] + self.newConfig.tau_Fu[1]) / 2,
#                                  (self.newConfig.tau_Fu[0] + self.newConfig.tau_Fu[1]) / 2]
#         self.newConfig.dtau_Fu = 1
#         self.newConfig.N = [(self.newConfig.N[0] + self.newConfig.N[1]) / 2,
#                             (self.newConfig.N[0] + self.newConfig.N[1]) / 2]
#         self.newConfig.dN = 1
#
#         self.pointsForTempOx = []
#         for currentTempOx in self.tempRangeOx:
#             self.newConfig.T_io = currentTempOx
#             unstablePointsMat = CCM.runPointAnalysis(self.newConfig)
#             for i in range(len(unstablePointsMat)):
#                 for j in range(len(unstablePointsMat[i])):
#                     for k in range(len(unstablePointsMat[i][j])):
#                         for l in range(len(unstablePointsMat[i][j][k])):
#                             for m in unstablePointsMat[i][j][k][l]:
#                                 tmpPoint = [currentTempOx, m]
#                                 self.pointsForTempOx.append(tmpPoint)
#
#         self.pointsForTempFu = []
#         for currentTempFu in self.tempRangeFu:
#             self.newConfig.T_if = currentTempFu
#             unstablePointsMat = CCM.runPointAnalysis(self.newConfig)
#             for i in range(len(unstablePointsMat)):
#                 for j in range(len(unstablePointsMat[i])):
#                     for k in range(len(unstablePointsMat[i][j])):
#                         for l in range(len(unstablePointsMat[i][j][k])):
#                             for m in unstablePointsMat[i][j][k][l]:
#                                 tmpPoint = [currentTempFu, m]
#                                 self.pointsForTempFu.append(tmpPoint)
#
#     def createTempVariationPlots(self):
#
#         workDir = os.getcwd()
#         path = workDir + "\\StabilityAnalysis"
#         if not (os.path.exists(path)):
#             os.mkdir(path)
#         path = workDir + "\\StabilityAnalysis\\" + self.newConfig.identifier
#         if not (os.path.exists(path)):
#             os.mkdir(path)
#
#         # Oxygen variation plot
#         plt.figure()
#         plt.grid()
#         plt.xlim(0, 1000)
#         plt.ylim(0, 5)
#         plt.xlabel("Frequency [Hz]")
#         plt.ylabel("Gain Factor")
#
#         x = [i[1][0] for i in self.pointsForTempOx]
#         y = [i[1][1] for i in self.pointsForTempOx]
#         col = [i[0] for i in self.pointsForTempOx]
#
#         sc = plt.scatter(x, y, c=col, marker="+")
#         cbar = plt.colorbar(sc)
#         cbar.set_label("Oxygen Temperature [K]")
#         plt.savefig(path + "\\OxygenTemperatureVariation")
#
#         # Fuel plot
#         plt.figure()
#         plt.grid()
#         plt.xlim(0, 1000)
#         plt.ylim(0, 5)
#         plt.xlabel("Frequency [Hz]")
#         plt.ylabel("Gain Factor")
#
#         x = [i[1][0] for i in self.pointsForTempFu]
#         y = [i[1][1] for i in self.pointsForTempFu]
#         col = [i[0] for i in self.pointsForTempFu]
#
#         sc = plt.scatter(x, y, c=col, marker="+")
#         cbar = plt.colorbar(sc)
#         cbar.set_label("Fuel Temperature [K]")
#         plt.savefig(path + "\\FuelTemperatureVariation")


######################################################################################################################################

# class geometricVariation:
#
#     def __init__(self, config, injector):
#         self.newConfig = config
#         self.init_injector = injector  # Copy of initial injector
#
#         # manifold changes the volume of the manifold
#         # lineLength changes length of the lines
#         # if both false -> changing the diameter of the lines
#         self.manifold = True
#         self.lineLength = False
#
#         self.run()
#         self.createGeometricVariationPlots()
#         self.createGeometricVariationBodePlots()
#
#     def run(self):
#
#         # Defining geometry variation range
#         startVariation = 0.75
#         endVariation = 2
#         self.variationVec = np.linspace(startVariation, endVariation, 10)
#
#         # Setting parameter values as average of input values
#         self.newConfig.tau_Ox = [(self.newConfig.tau_Ox[0] + self.newConfig.tau_Ox[1]) / 2,
#                                  (self.newConfig.tau_Ox[0] + self.newConfig.tau_Ox[1]) / 2]
#         self.newConfig.dtau_Ox = 1
#         self.newConfig.tau_Fu = [(self.newConfig.tau_Fu[0] + self.newConfig.tau_Fu[1]) / 2,
#                                  (self.newConfig.tau_Fu[0] + self.newConfig.tau_Fu[1]) / 2]
#         self.newConfig.dtau_Fu = 1
#         self.newConfig.N = [(self.newConfig.N[0] + self.newConfig.N[1]) / 2,
#                             (self.newConfig.N[0] + self.newConfig.N[1]) / 2]
#         self.newConfig.dN = 1
#
#         self.pointsForTempOx = []
#         self.solVectorsOx = []
#         for currentVarOx in self.variationVec:
#             print(currentVarOx)
#             newInjector = self.writeChangedGeometrieToExcel("Oxygen", currentVarOx, self.lineLength, self.manifold)
#
#             unstablePointsMat, solVec = CCM.runPointAnalysis(self.newConfig, newInjector)
#
#             self.solVectorsOx.append(solVec)
#             for i in range(len(unstablePointsMat)):
#                 for j in range(len(unstablePointsMat[i])):
#                     for k in range(len(unstablePointsMat[i][j])):
#                         for l in range(len(unstablePointsMat[i][j][k])):
#                             for m in unstablePointsMat[i][j][k][l]:
#                                 tmpPoint = [currentVarOx, m]
#                                 self.pointsForTempOx.append(tmpPoint)
#
#         # startVariation = 0.9
#         # endVariation = 2
#         # self.variationVec = np.linspace(startVariation, endVariation, 10)
#
#         self.pointsForTempFu = []
#         self.solVectorsFu = []
#         for currentVarFu in self.variationVec:
#             print(currentVarFu)
#             newInjector = self.writeChangedGeometrieToExcel("Fuel", currentVarFu, self.lineLength, self.manifold)
#             unstablePointsMat, solVec = CCM.runPointAnalysis(self.newConfig, newInjector)
#             self.solVectorsFu.append(solVec)
#             for i in range(len(unstablePointsMat)):
#                 for j in range(len(unstablePointsMat[i])):
#                     for k in range(len(unstablePointsMat[i][j])):
#                         for l in range(len(unstablePointsMat[i][j][k])):
#                             for m in unstablePointsMat[i][j][k][l]:
#                                 tmpPoint = [currentVarFu, m]
#                                 self.pointsForTempFu.append(tmpPoint)
#
#     def createGeometricVariationPlots(self):
#
#         workDir = os.getcwd()
#         path = workDir + "\\StabilityAnalysis"
#         if not (os.path.exists(path)):
#             os.mkdir(path)
#         path = workDir + "\\StabilityAnalysis\\" + self.newConfig.identifier
#         if not (os.path.exists(path)):
#             os.mkdir(path)
#
#         # Oxygen variation plot
#         plt.figure()
#         plt.grid()
#         plt.xlim(0, 1000)
#         plt.ylim(0, 5)
#         plt.xlabel("Frequency [Hz]")
#         plt.ylabel("Gain Factor")
#
#         x = [i[1][0] for i in self.pointsForTempOx]
#         y = [i[1][1] for i in self.pointsForTempOx]
#         col = [i[0] for i in self.pointsForTempOx]
#
#         if (self.manifold):
#             element = "Manifold_Volume_"
#         else:
#             if (self.lineLength):
#                 element = "Line_Length_"
#             else:
#                 element = "Line_Diameter_"
#
#         sc = plt.scatter(x, y, c=col, marker="+")
#         cbar = plt.colorbar(sc)
#         cbar.set_label(element + "Variation")
#         plt.savefig(path + "\\Oxygen_" + element + "Variation")
#
#         # Fuel plot
#         plt.figure()
#         plt.grid()
#         plt.xlim(0, 1000)
#         plt.ylim(0, 5)
#         plt.xlabel("Frequency [Hz]")
#         plt.ylabel("Gain Factor")
#
#         x = [i[1][0] for i in self.pointsForTempFu]
#         y = [i[1][1] for i in self.pointsForTempFu]
#         col = [i[0] for i in self.pointsForTempFu]
#
#         sc = plt.scatter(x, y, c=col, marker="+")
#         cbar = plt.colorbar(sc)
#         cbar.set_label(element + "Variation")
#         plt.savefig(path + "\\Fuel_" + element + "Variation")
#
#     def createGeometricVariationBodePlots(self):
#
#         workDir = os.getcwd()
#         path = workDir + "\\StabilityAnalysis"
#         if not (os.path.exists(path)):
#             os.mkdir(path)
#         path = workDir + "\\StabilityAnalysis\\" + self.newConfig.identifier
#         if not (os.path.exists(path)):
#             os.mkdir(path)
#
#         if (self.manifold):
#             element = "Manifold_Volume_"
#         else:
#             if (self.lineLength):
#                 element = "Line_Length_"
#             else:
#                 element = "Line_Diameter_"
#
#         df = self.newConfig.df
#         fmin = self.newConfig.freq_range[0]
#         fmax = self.newConfig.freq_range[1]
#         freqs = np.arange(fmin, fmax + df * 0.1, df)
#
#         # Oxygen variation plot
#         plt.figure()
#         sb1 = plt.subplot(2, 1, 1)
#         sb1.grid()
#         sb1.set_ylabel("Gain")
#
#         sb2 = plt.subplot(2, 1, 2)
#         sb2.grid()
#         sb2.set_xlabel("Frequency [Hz]")
#         sb2.set_ylabel("Phase")
#
#         for j in range(len(self.solVectorsOx)):
#             sb1.plot(freqs, abs(self.solVectorsOx[j]))
#             sb1.set_xlim(self.newConfig.freq_range[0], self.newConfig.freq_range[1])
#             yLim = [abs(i) for i in self.solVectorsOx[j] if not np.isnan(i)]
#             sb1.set_ylim(0, 5)
#
#             phi = np.arctan2(np.imag(self.solVectorsOx[j]), np.real(self.solVectorsOx[j]))
#             phi = phi * 180 / np.pi
#             sb2.plot(freqs, phi)
#             sb2.set_xlim(self.newConfig.freq_range[0], self.newConfig.freq_range[1])
#
#         plt.savefig(path + "\\Bode_Diagram_Oxygen_" + element + "Variation")
#
#         # Fuel variation plot
#         plt.figure()
#         sb1 = plt.subplot(2, 1, 1)
#         sb1.grid()
#         sb1.set_ylabel("Gain")
#
#         sb2 = plt.subplot(2, 1, 2)
#         sb2.grid()
#         sb2.set_xlabel("Frequency [Hz]")
#         sb2.set_ylabel("Phase")
#
#         for G in self.solVectorsFu:
#             sc = sb1.plot(freqs, abs(G))
#             sb1.set_xlim(self.newConfig.freq_range[0], self.newConfig.freq_range[1])
#             yLim = [abs(i) for i in G if not np.isnan(i)]
#             sb1.set_ylim(0, max(yLim) * 1.1)
#
#             phi = np.arctan2(np.imag(G), np.real(G))
#             phi = phi * 180 / np.pi
#             sb2.plot(freqs, phi)
#             sb2.set_xlim(self.newConfig.freq_range[0], self.newConfig.freq_range[1])
#
#         plt.savefig(path + "\\Bode_Diagram_Fuel_" + element + "Variation")
#
#     def writeChangedGeometrieToExcel(self, fluid, change, length, manifold):
#
#         loc = ("Geometry_Injection_Configurations" + "\\" + self.newConfig.injector_configuration + ".xlsx")
#         wb = xlrd.open_workbook(loc)
#         injector_file = wb.sheet_by_index(0)
#
#         # Calculation for the Casiano Injector mode has to be delayed
#         casMode = False
#         if self.newConfig.casMode:
#             self.newConfig.casMode = False
#             casMode = True
#
#         newInjector = injector_configuration(injector_file, self.newConfig)
#
#         if (fluid == "Oxygen"):
#             for i in range(len(newInjector.LOX.elements)):
#                 if (manifold):
#                     if (newInjector.LOX.elements[i].type == "Manifold"):
#                         newInjector.LOX.elements[i].volume *= change
#                 else:
#                     if (newInjector.LOX.elements[i].type == "Line"):
#                         if (length):
#                             newInjector.LOX.elements[i].length *= change
#                         else:
#                             newInjector.LOX.elements[i].diameter *= change
#                             newInjector.LOX.elements[i].area *= change ** 2
#
#         elif (fluid == "Fuel"):
#             for i in range(len(newInjector.Fuel.elements)):
#                 if (manifold):
#                     if (newInjector.Fuel.elements[i].type == "Manifold"):
#                         newInjector.Fuel.elements[i].volume *= change
#                 else:
#                     if (newInjector.Fuel.elements[i].type == "Line"):
#                         if (length):
#                             newInjector.Fuel.elements[i].length *= change
#                         else:
#                             newInjector.Fuel.elements[i].diameter *= change
#                             newInjector.Fuel.elements[i].area *= change ** 2
#         else:
#             print("ERROR")
#             sys.exit()
#
#         newInjector.assignFluidTypes()
#         newInjector.assignTemperaturesAndMassflowsToElements()
#         newInjector.calculatePressureInElements()
#         newInjector.calculateAcousticParameters()
#
#         if (casMode):
#             self.newConfig.casMode = True
#             newInjector.reduceInjectorHeadToSingleInjector()
#
#         return newInjector


# class pressureDropIteration:
#
#     def __init__(self, config):
#         self.newConfig = config
#         self.newConfig.freq_range[0] = 200
#         self.newConfig.freq_range[1] = 1000
#         self.opPressOx = self.newConfig.p_io
#         self.opPressFu = self.newConfig.p_if
#         self.changeDiameterForMassFlow = False
#         self.linePosOx = 1
#         self.linePosFu = 3
#         self.upperThresholdOxygenPDropInit = config.p_cc * 1.5
#         self.lowerThresholdOxygenPDropInit = config.p_cc * 1.03
#         self.fuelPressureDrops = np.linspace(config.p_cc * 1.03, config.p_cc * 2.5, 25)
#
#         self.run()
#         self.createPressureDropStabilityPlot()
#         self.saveResultToFile()
#
#     def checkStability(self, oxyPres, fuelPres):
#
#         oldDeltaPOx = self.opPressOx - self.newConfig.p_cc
#         newDeltaPOx = oxyPres - self.newConfig.p_cc
#
#         if (self.changeDiameterForMassFlow):
#             try:
#                 newInjector = self.changedGeometrieForConstMassFlow(oldDeltaPOx, newDeltaPOx)
#                 self.newConfig.p_io = oxyPres
#                 self.newConfig.p_if = fuelPres
#                 unstablePointsMat = CCM.runPointAnalysis(self.newConfig, newInjector)
#             except:
#                 return False
#         else:
#             self.newConfig.p_io = oxyPres
#             self.newConfig.p_if = fuelPres
#
#             try:
#                 unstablePointsMat = CCM.runPointAnalysis(self.newConfig)
#             except:
#                 print("Stability Check Failure")
#                 return False
#
#         print("Stability Check Successful")
#         stable = True
#         for i in range(len(unstablePointsMat)):
#             for j in range(len(unstablePointsMat[i])):
#                 for k in range(len(unstablePointsMat[i][j])):
#                     for l in range(len(unstablePointsMat[i][j][k])):
#                         if not (len(unstablePointsMat[i][j][k][l]) == 0):
#                             stable = False
#
#         return stable
#
#     def getPressureDropInPercent(self, pressure):
#         return ((pressure - self.newConfig.p_cc) / self.newConfig.p_cc * 100)
#
#     def run(self):
#
#         self.stabilityMap = []
#         for currentFuelPDrop in self.fuelPressureDrops:
#             print("Fuel Pressure Frop: " + str(self.getPressureDropInPercent(currentFuelPDrop)))
#             upperBound = self.upperThresholdOxygenPDropInit
#             lowerBound = self.lowerThresholdOxygenPDropInit
#
#             upStable = self.checkStability(upperBound, currentFuelPDrop)
#             lowStable = self.checkStability(lowerBound, currentFuelPDrop)
#
#             if (
#                     upStable and not lowStable):  # Upper Limit stable -> stability iteration possible,  low Limit stable -> always stable, no iteration necessary
#                 midPressureDrop = (upperBound + lowerBound) / 2
#                 epsilon = upperBound - midPressureDrop
#
#                 while (epsilon > 0.01):
#                     midStable = self.checkStability(midPressureDrop, currentFuelPDrop)
#
#                     if (midStable):
#                         upperBound = midPressureDrop
#                         midPressureDrop = (upperBound + lowerBound) / 2
#                         epsilon = upperBound - midPressureDrop
#                     else:
#                         lowerBound = midPressureDrop
#                         midPressureDrop = (upperBound + lowerBound) / 2
#                         epsilon = upperBound - midPressureDrop
#
#                 tmpObj = [self.getPressureDropInPercent(upperBound), self.getPressureDropInPercent(currentFuelPDrop)]
#                 self.stabilityMap.append(tmpObj)
#
#             elif (upStable and lowStable):
#                 print("Inherently stable")
#                 tmpObj = [5, self.getPressureDropInPercent(currentFuelPDrop)]
#                 self.stabilityMap.append(tmpObj)
#
#             else:
#                 print("Inherently unstable")
#                 tmpObj = [np.NaN, self.getPressureDropInPercent(currentFuelPDrop)]
#                 self.stabilityMap.append(tmpObj)
#
#     def createPressureDropStabilityPlot(self):
#
#         plt.figure(figsize=(4.5, 4), dpi=250)
#         x = [i[0] for i in self.stabilityMap]
#         y = [i[1] for i in self.stabilityMap]
#         plt.plot(x, y, c="r")
#         x = self.getPressureDropInPercent(self.opPressOx)
#         y = self.getPressureDropInPercent(self.opPressFu)
#         plt.plot(x, y, c="b", marker="x")
#         plt.xlabel("Oxygen pressure drop $\Delta P_{Ox}$ in %")
#         plt.ylabel("Fuel pressure drop $\Delta P_{Fu}$ in %")
#         plt.xlim(0, self.getPressureDropInPercent(self.upperThresholdOxygenPDropInit))
#         plt.ylim(0, self.getPressureDropInPercent(self.fuelPressureDrops[-1]))
#         plt.show()
#         plt.savefig(
#             os.getcwd() + "\\StabilityAnalysis\\" + self.newConfig.identifier + "\\PressureDropStabilityPlot.png")
#
#     def saveResultToFile(self):
#         workDir = os.getcwd()
#         path = workDir + "\\StabilityAnalysis"
#         if not (os.path.exists(path)):
#             os.mkdir(path)
#         path = path + "\\" + self.newConfig.identifier
#         if not (os.path.exists(path)):
#             os.mkdir(path)
#         filePath = path + "\\PressureDropIteration.dat"
#
#         # Opening and writing header for the resultFile
#         if (os.path.exists(filePath)):
#             fileTxt = open(filePath, "a")
#         else:
#             try:
#                 fileTxt = open(filePath, "wt")
#                 fileTxt.write("DeltaP_Ox\tDeltaP_Fu\n")
#             except:
#                 fileTxt = open(filePath, "w")
#
#         for i in self.stabilityMap:
#             fileTxt.write(str(format(np.round(i[0], 5), ".5f")) + "\t" + str(format(np.round(i[1], 5), ".5f")) + "\n")
#
#         fileTxt.close()
#
#     def changedGeometrieForConstMassFlow(self, oldDeltaP, newDeltaP):
#
#         loc = ("Geometry_Injection_Configurations" + "\\" + self.newConfig.injector_configuration + ".xlsx")
#         wb = xlrd.open_workbook(loc)
#         injector_file = wb.sheet_by_index(0)
#
#         # Calculation for the Casiano Injector mode has to be delayed
#         casMode = False
#         if (self.newConfig.casMode):
#             self.newConfig.casMode = False
#             casMode = True
#
#         print("Test1")
#
#         self.newConfig.p_io = self.opPressOx
#         self.newConfig.p_if = self.opPressFu
#
#         newInjector = injector_configuration(injector_file, self.newConfig)
#
#         massFlowSinOx = self.newConfig.mDot_O2 / newInjector.LOX.elements[-1].numOfElements
#         massFlowSinFu = self.newConfig.mDot_F / newInjector.Fuel.elements[-1].numOfElements
#
#         counter = 0
#         for i in range(len(newInjector.LOX.elements)):
#             if (newInjector.LOX.elements[i].type == "Line"):
#                 counter += 1
#                 if (counter == self.linePosOx):
#                     positionOx = i
#                     break
#         print("Test2")
#
#         change = np.sqrt(massFlowSinOx / (newInjector.LOX.elements[positionOx].area
#                                           * np.sqrt(
#                     2 * newInjector.LOX.elements[positionOx].rho * newInjector.LOX.elements[
#                         positionOx - 1].dp * 10 ** 5 * newDeltaP / oldDeltaP)))
#
#         print(change)
#
#         counter = 0
#         for i in range(len(newInjector.LOX.elements)):
#             if (newInjector.LOX.elements[i].type == "Line"):
#                 counter += 1
#                 if (counter == self.linePosOx):
#                     newInjector.LOX.elements[i].diameter *= change
#                     newInjector.LOX.elements[i].area *= change ** 2
#                     break
#
#         # counter = 0
#         # for i in range(len(newInjector.Fuel.elements)):
#         #     if (newInjector.Fuel.elements[i].type == "Line"):
#         #         counter += 1
#         #         if (counter == self.linePosOx):
#         #             newInjector.Fuel.elements[i].diameter *= change
#         #             newInjector.Fuel.elements[i].area *= change**2
#         #             break
#
#         newInjector.assignFluidTypes()
#         newInjector.assignTemperaturesAndMassflowsToElements()
#         newInjector.calculatePressureInElements()
#         newInjector.calculateAcousticParameters()
#
#         if (casMode):
#             self.newConfig.casMode = True
#             newInjector.reduceInjectorHeadToSingleInjector()
#
#         return newInjector
