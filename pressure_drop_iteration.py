""" Stability research mode - pressure drop iteration

Created on Wed Jun 9 07:36:02 2021

"""

import numpy as np
from point_analysis import run_point_analysis
from tqdm import tqdm
import matplotlib.pyplot as plt
import os
import pathlib
import csv


RESULT_DIRECTORY = 'StabilityFieldPrediction'


class PressureDropIteration:
    """ Pressure drop iteration, or stability field prediction.

    See the master thesis of Geuking, DLR-LA-RAK-HF-RP-110.
    """

    def __init__(self, config, change_diameter_for_mass_flow=False):
        """ Creates an object of :class:`PressureDropIteration` and starts the calculation, then creates the
        corresponding plot and stores the results in a file.

        :param config: Containing all initial settings
        :param change_diameter_for_mass_flow: Determines of the diameter for the mass flow is changed, default False

        :type config: :class:`Config.ParseConfig`
        :type change_diameter_for_mass_flow: bool
        :rtype: None
        """

        # Creating a copy of the initial config
        self.newConfig = config

        # Re-defining the frequency range manually
        self.newConfig.freq_range = np.array([200, 1000])

        self.opPressOx = self.newConfig.p_io
        self.opPressFu = self.newConfig.p_if
        self.changeDiameterForMassFlow = change_diameter_for_mass_flow
        self.linePosOx = 1  # ?
        self.linePosFu = 3  # ?
        self.upperThresholdOxygenPDropInit = config.p_cc * 1.5  # ?
        self.lowerThresholdOxygenPDropInit = config.p_cc * 1.03  # ?
        self.fuelPressureDrops = np.linspace(config.p_cc * 1.03, config.p_cc * 2.5, 25)  # ?

        # Create an object to monitor the progress. This is just for monitoring and display purposes
        monitor_progress = tqdm(total=len(self.fuelPressureDrops),
                                desc='Pressure drop iteration for {}'.format(self.newConfig.identifier),
                                leave=True, bar_format='{desc:<65}| {bar:25}| {n_fmt}/{total_fmt}')

        # Empty stability map to collect and store results
        self.stabilityMap = []

        # Calculate pressure drop iteration
        for currentFuelPDrop in self.fuelPressureDrops:

            # For a current fuel pressure drop, find (if possible) the corresponding oxygen pressure drop
            # where the system is stable

            # print("Fuel Pressure Frop: " + str(self.get_pressure_drop_in_percent(currentFuelPDrop)))

            # Defining the current boundaries (for the oxygen pressure drop)
            upperBound = self.upperThresholdOxygenPDropInit
            lowerBound = self.lowerThresholdOxygenPDropInit

            # Check if the upper and if the lower boundary result in a stable behavior
            upStable = self.check_stability(upperBound, currentFuelPDrop)
            lowStable = self.check_stability(lowerBound, currentFuelPDrop)

            if (upStable and not lowStable):
                # Upper Limit stable -> stability iteration possible,
                # low Limit stable -> always stable, no iteration necessary
                midPressureDrop = (upperBound + lowerBound) / 2
                epsilon = upperBound - midPressureDrop

                while epsilon > 0.01:
                    midStable = self.check_stability(midPressureDrop, currentFuelPDrop)

                    if midStable:
                        upperBound = midPressureDrop
                        midPressureDrop = (upperBound + lowerBound) / 2
                        epsilon = upperBound - midPressureDrop
                    else:
                        lowerBound = midPressureDrop
                        midPressureDrop = (upperBound + lowerBound) / 2
                        epsilon = upperBound - midPressureDrop

                tmpObj = [self.get_pressure_drop_in_percent(upperBound), self.get_pressure_drop_in_percent(currentFuelPDrop)]
                self.stabilityMap.append(tmpObj)

            elif upStable and lowStable:
                # print("Inherently stable")
                # TODO: Why '5'? -> 5 is the oxygen pressure drop in percent
                tmpObj = [5, self.get_pressure_drop_in_percent(currentFuelPDrop)]
                self.stabilityMap.append(tmpObj)

            else:
                # This equals 'not upStable and lowStable'
                # print("Inherently unstable")
                tmpObj = [np.NaN, self.get_pressure_drop_in_percent(currentFuelPDrop)]
                self.stabilityMap.append(tmpObj)

            # Debug
            # print(tmpObj)

            # Update monitoring, just for display purposes
            monitor_progress.update()

        # Close monitoring, just for display purposes
        monitor_progress.close()

        # Create path to store results
        path = os.path.join(self.newConfig.result_path, RESULT_DIRECTORY)
        pathlib.Path(path).mkdir(parents=True)

        # Create plot
        plt.figure(figsize=(4.5, 4), dpi=250)
        x = [i[0] for i in self.stabilityMap]
        y = [i[1] for i in self.stabilityMap]
        plt.plot(x, y, c="r")
        x = self.get_pressure_drop_in_percent(self.opPressOx)
        y = self.get_pressure_drop_in_percent(self.opPressFu)
        plt.plot(x, y, c="b", marker="x")
        plt.xlabel("Oxygen pressure drop $\Delta P_{Ox}$ in %")
        plt.ylabel("Fuel pressure drop $\Delta P_{Fu}$ in %")
        plt.xlim(0, self.get_pressure_drop_in_percent(self.upperThresholdOxygenPDropInit))
        plt.ylim(0, self.get_pressure_drop_in_percent(self.fuelPressureDrops[-1]))
        plt.show()
        plt.savefig(os.path.join(path, 'PressureDropStabilityPlot.png'))

        # Store results in file

        # workDir = os.getcwd()
        # path = workDir + "\\StabilityAnalysis"
        # if not (os.path.exists(path)):
        #     os.mkdir(path)
        # path = path + "\\" + self.newConfig.identifier
        # if not (os.path.exists(path)):
        #     os.mkdir(path)
        # filePath = path + "\\PressureDropIteration.dat"
        #
        # Opening and writing header for the resultFile
        # if (os.path.exists(filePath)):
        #     fileTxt = open(filePath, "a")
        # else:
        #     try:
        #         fileTxt = open(filePath, "wt")
        #         fileTxt.write("DeltaP_Ox\tDeltaP_Fu\n")
        #     except:
        #         fileTxt = open(filePath, "w")
        #
        # for i in self.stabilityMap:
        #     fileTxt.write(str(format(np.round(i[0], 5), ".5f")) + "\t" + str(format(np.round(i[1], 5), ".5f")) + "\n")
        #
        # fileTxt.close()

        with open(os.path.join(path, 'PressureDropIteration.dat'), 'w+', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=['Delta_P_Ox', 'Delta_P_Fu'], delimiter='\t')
            writer.writeheader()
            for delta_p_ox, delta_p_fu in self.stabilityMap:
                writer.writerow({'Delta_P_Ox': np.round(delta_p_ox, 5), 'Delta_P_Fu': np.round(delta_p_fu, 5)})

    def check_stability(self, oxyPres, fuelPres):
        """ Determines if a given combination of an oxygen pressure drop and a fuel pressure drop is stable.

        :param oxyPres: Oxygen pressure drop
        :param fuelPres: Fuel pressure drop
        :return: Stable, or unstable

        :type oxyPres: float
        :type fuelPres: float
        :rtype: bool
        """

        oldDeltaPOx = self.opPressOx - self.newConfig.p_cc
        newDeltaPOx = oxyPres - self.newConfig.p_cc

        if self.changeDiameterForMassFlow:
            try:
                newInjector = self.changed_geometry_for_constant_mass_flow(oldDeltaPOx, newDeltaPOx)
                self.newConfig.p_io = oxyPres
                self.newConfig.p_if = fuelPres
                unstablePointsMat, _ = run_point_analysis(self.newConfig, newInjector,
                                                          self.newConfig.combustion_chamber)
            except:
                return False

        else:
            # Changing the injector properties, which concerns the injector - the injector has to be reloaded
            self.newConfig.p_io = oxyPres
            self.newConfig.p_if = fuelPres

            try:
                # Using the initial injector but reload because the pressure was changed before
                unstablePointsMat, _ = run_point_analysis(self.newConfig, self.newConfig.reload_injector(),
                                                       self.newConfig.combustion_chamber)
            except:
                # print("Stability Check Failure")
                return False

        # print("Stability Check Successful")
        stable = True
        for i in range(len(unstablePointsMat)):
            for j in range(len(unstablePointsMat[i])):
                for k in range(len(unstablePointsMat[i][j])):
                    for l in range(len(unstablePointsMat[i][j][k])):
                        if not (len(unstablePointsMat[i][j][k][l]) == 0):
                            stable = False

        return stable

    def get_pressure_drop_in_percent(self, pressure):
        """ Calculates the pressure drop in percent between a given pressure and the combustion chamber pressure.

        :param pressure: Current pressure
        :return: Pressure drop in percent

        :type pressure: float
        :rtype: float
        """
        return (pressure - self.newConfig.p_cc) / self.newConfig.p_cc * 100

    def changed_geometry_for_constant_mass_flow(self, oldDeltaP, newDeltaP):
        """

        :param oldDeltaP:
        :param newDeltaP:
        :return:
        """

        # loc = ("Geometry_Injection_Configurations" + "\\" + self.newConfig.injector_configuration + ".xlsx")
        # wb = xlrd.open_workbook(loc)
        # injector_file = wb.sheet_by_index(0)

        # TODO: ? (Confusing logic, see end of function) -> just use self.newConfig.casMode, which is used in
        #  injector.recalculate_injector() anyway.
        # Calculation for the Casiano Injector mode has to be delayed
        # casMode = False
        # if self.newConfig.casMode:
        #     self.newConfig.casMode = False
        #     casMode = True

        # print("Test1")

        self.newConfig.p_io = self.opPressOx
        self.newConfig.p_if = self.opPressFu

        # newInjector = injector_configuration(injector_file, self.newConfig)
        newInjector = self.newConfig.reload_injector()

        massFlowSinOx = self.newConfig.mDot_O2 / newInjector.LOX.elements[-1].numOfElements
        massFlowSinFu = self.newConfig.mDot_F / newInjector.Fuel.elements[-1].numOfElements

        counter = 0
        for i in range(len(newInjector.LOX.elements)):
            if newInjector.LOX.elements[i].type == "Line":
                counter += 1
                if counter == self.linePosOx:
                    positionOx = i
                    break
        # print("Test2")

        change = np.sqrt(massFlowSinOx / (newInjector.LOX.elements[positionOx].area *
                                          np.sqrt(2 * newInjector.LOX.elements[positionOx].rho *
                                                  newInjector.LOX.elements[positionOx - 1].dp * 10 ** 5 * newDeltaP /
                                                  oldDeltaP)))

        # print(change)

        counter = 0
        for i in range(len(newInjector.LOX.elements)):
            if (newInjector.LOX.elements[i].type == "Line"):
                counter += 1
                if (counter == self.linePosOx):
                    newInjector.LOX.elements[i].diameter *= change
                    newInjector.LOX.elements[i].area *= change ** 2
                    break

        # Note: This was already commented earlier, has no function.
        # counter = 0
        # for i in range(len(newInjector.Fuel.elements)):
        #     if (newInjector.Fuel.elements[i].type == "Line"):
        #         counter += 1
        #         if (counter == self.linePosOx):
        #             newInjector.Fuel.elements[i].diameter *= change
        #             newInjector.Fuel.elements[i].area *= change**2
        #             break

        # newInjector.assignFluidTypes()
        # newInjector.assignTemperaturesAndMassflowsToElements()
        # newInjector.calculatePressureInElements()
        # newInjector.calculateAcousticParameters()
        #
        # if (casMode):
        #     self.newConfig.casMode = True
        #     newInjector.reduceInjectorHeadToSingleInjector()

        # Recalculating properties for the changed injector geometry
        newInjector.recalculate_injector()

        return newInjector
