""" Stability research mode - injector geometry variation

Created on Mon Jun 7 14:32:21 2021

"""

import numpy as np
import pathlib
import os
import matplotlib.pyplot as plt
from point_analysis import run_point_analysis
from tqdm import tqdm

RESULT_DIRECTORY = 'InjectorGeometryVariation'


def injector_geometry_variation(config, range_factor_low, range_factor_high, steps, vary_line_diameter, vary_line_length,
                                vary_manifold):
    """ Varies the injector diameter, the injector line length, and/or the injector manifold volume.

    :param vary_line_length: Determines if the line length should be varied
    :param vary_line_diameter: Determined if the line diameter should be varied
    :param vary_manifold: Determines if the manifold volume should be varied
    :param config: A :class:`Config.ParseConfig` object containing the initial settings
    :param range_factor_low: Lower boundary factor for the variation
    :param range_factor_high: Upper boundary factor for the variation
    :param steps: Amount of variation steps between the lower and upper boundary factors

    :type config: :class:`Config.ParseConfig`
    :type range_factor_low: float
    :type range_factor_high: float
    :type steps: int
    :type vary_line_diameter: bool
    :type vary_line_length: bool
    :type vary_manifold: bool
    :rtype: None
    """

    new_config = config
    # Copy of the initial, default injector
    default_injector = config.injector

    # Variation spectrum; conversion to int because np.linspace is expecting an integer
    variationVec = np.linspace(range_factor_low, range_factor_high, int(steps))

    # Because the stability analysis should not be evaluated for multiple values of tau_Ox, tau_Fu and N, new values are
    # set to be the average of the input values
    new_config.tau_Ox = [(min(new_config.tau_Ox) + max(new_config.tau_Ox)) / 2, (min(new_config.tau_Ox) + max(new_config.tau_Ox)) / 2]
    new_config.dtau_Ox = 1
    new_config.tau_Fu = [(min(new_config.tau_Fu) + max(new_config.tau_Fu)) / 2, (min(new_config.tau_Fu) + max(new_config.tau_Fu)) / 2]
    new_config.dtau_Fu = 1
    new_config.N = [(min(new_config.N) + max(new_config.N)) / 2, (min(new_config.N) + max(new_config.N)) / 2]
    new_config.dN = 1

    # Just for monitoring and display purposed, creates a progress bar
    monitor_progress = tqdm(total=steps*2, desc='Injector geometry variation for {}'.format(new_config.identifier),
                            leave=True, bar_format='{desc:<65}| {bar:25}| {n_fmt}/{total_fmt}')

    # Calculate the solution for the changing geometry. This is done for 'Oxygen' and 'Fuel' separately.
    # For each value in variationVec, a new injector is created and the solution is calculated. Based on all collected
    # solutions, the plots are created afterwards
    # Oxygen
    pointsForTempOx = []
    solVectorsOx = []
    for currentVarOx in variationVec:
        newInjector = change_injector_geometry(config, "Oxygen", currentVarOx,
                                               vary_line_length, vary_line_diameter, vary_manifold)

        unstablePointsMat, solVec = run_point_analysis(new_config, newInjector, new_config.combustion_chamber)
        solVectorsOx.append(solVec)

        for i in range(len(unstablePointsMat)):
            for j in range(len(unstablePointsMat[i])):
                for k in range(len(unstablePointsMat[i][j])):
                    for l in range(len(unstablePointsMat[i][j][k])):
                        for m in unstablePointsMat[i][j][k][l]:
                            tmpPoint = [currentVarOx, m]
                            pointsForTempOx.append(tmpPoint)
            # Updates the progress bar, only for monitoring purposes
            monitor_progress.update()

    # Fuel
    pointsForTempFu = []
    solVectorsFu = []
    for currentVarFu in variationVec:
        newInjector = change_injector_geometry(config, "Fuel", currentVarFu,
                                               vary_line_length, vary_line_diameter, vary_manifold)
        unstablePointsMat, solVec = run_point_analysis(new_config, newInjector, new_config.combustion_chamber)
        solVectorsFu.append(solVec)

        for i in range(len(unstablePointsMat)):
            for j in range(len(unstablePointsMat[i])):
                for k in range(len(unstablePointsMat[i][j])):
                    for l in range(len(unstablePointsMat[i][j][k])):
                        for m in unstablePointsMat[i][j][k][l]:
                            tmpPoint = [currentVarFu, m]
                            pointsForTempFu.append(tmpPoint)
            # Updates the progress bar, only for monitoring purposes
            monitor_progress.update()

    # Finish monitoring, just for display purposes
    monitor_progress.close()

    # Create plots and store results
    # Create the directory to store results
    path = os.path.join(new_config.result_path, RESULT_DIRECTORY)
    pathlib.Path(path).mkdir(parents=True)

    # Determine which element(s) were changed, for naming the plots correctly
    element = ''
    if vary_manifold:
        element += "Manifold_Volume_"
    if vary_line_length:
        element += "Line_Length_"
    if vary_line_diameter:
        element += "Line_Diameter_"

    # Oxygen variation plot
    plt.figure()
    plt.grid()
    # plt.xlim(0, 1000)
    plt.xlim(min(new_config.freq_range), max(new_config.freq_range))
    plt.ylim(0, 5)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Gain Factor")
    x = [i[1][0] for i in pointsForTempOx]
    y = [i[1][1] for i in pointsForTempOx]
    col = [i[0] for i in pointsForTempOx]
    sc = plt.scatter(x, y, c=col, marker="+")
    cbar = plt.colorbar(sc)
    cbar.set_label(element + "Variation")
    # plt.savefig(path + "\\Oxygen_" + element + "Variation")
    plt.savefig(os.path.join(path, 'Oxygen_{}_Variation'.format(element)))

    # Fuel variation plot
    plt.figure()
    plt.grid()
    # plt.xlim(0, 1000)
    plt.xlim(min(new_config.freq_range), max(new_config.freq_range))
    plt.ylim(0, 5)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Gain Factor")
    x = [i[1][0] for i in pointsForTempFu]
    y = [i[1][1] for i in pointsForTempFu]
    col = [i[0] for i in pointsForTempFu]
    sc = plt.scatter(x, y, c=col, marker="+")
    cbar = plt.colorbar(sc)
    cbar.set_label(element + "Variation")
    # plt.savefig(path + "\\Fuel_" + element + "Variation")
    plt.savefig(os.path.join(path, 'Fuel_{}_Variation'.format(element)))

    # Create bode plots
    frequency_range = np.arange(min(new_config.freq_range), max(new_config.freq_range) + 0.1 * new_config.df, new_config.df)

    # Oxygen bode plot
    plt.figure()
    sb1 = plt.subplot(2, 1, 1)
    sb1.grid()
    sb1.set_ylabel("Gain")
    sb2 = plt.subplot(2, 1, 2)
    sb2.grid()
    sb2.set_xlabel("Frequency [Hz]")
    sb2.set_ylabel("Phase")
    for j in range(len(solVectorsOx)):
        sb1.plot(frequency_range, abs(solVectorsOx[j]))
        sb1.set_xlim(min(new_config.freq_range), max(new_config.freq_range))
        yLim = [abs(i) for i in solVectorsOx[j] if not np.isnan(i)]
        sb1.set_ylim(0, 5)
        phi = np.arctan2(np.imag(solVectorsOx[j]), np.real(solVectorsOx[j]))
        phi = phi * 180 / np.pi
        sb2.plot(frequency_range, phi)
        sb2.set_xlim(new_config.freq_range[0], new_config.freq_range[1])
    plt.savefig(os.path.join(path, 'Bode_Diagram_Oxygen_{}_Variatio'.format(element)))

    # Fuel bode plot
    plt.figure()
    sb1 = plt.subplot(2, 1, 1)
    sb1.grid()
    sb1.set_ylabel("Gain")
    sb2 = plt.subplot(2, 1, 2)
    sb2.grid()
    sb2.set_xlabel("Frequency [Hz]")
    sb2.set_ylabel("Phase")
    for G in solVectorsFu:
        sc = sb1.plot(frequency_range, abs(G))
        sb1.set_xlim(min(new_config.freq_range), max(new_config.freq_range))
        yLim = [abs(i) for i in G if not np.isnan(i)]
        sb1.set_ylim(0, max(yLim) * 1.1)
        phi = np.arctan2(np.imag(G), np.real(G))
        phi = phi * 180 / np.pi
        sb2.plot(frequency_range, phi)
        sb2.set_xlim(new_config.freq_range[0], new_config.freq_range[1])
    plt.savefig(os.path.join(path, 'Bode_Diagram_Fuel_{}_Variation'.format(element)))


def change_injector_geometry(config, fluid, change, line_length, line_diameter, manifold):
    """ Creates a injector with the changed geometry.

    :param config: Current configuration
    :param fluid: Defines the fluid in the injector, 'Oxygen' or 'Fuel'
    :param change: Current change factor
    :param line_length: Determines if the line length should be varied
    :param line_diameter: Determines if the line diameter should be varied
    :param manifold: Determined if te manifold volume should be varied
    :return: New injector with a changed geometry

    :type config: :class:`Config.ParseConfig`
    :type fluid: str
    :type change: float
    :type line_length: bool
    :type line_diameter: bool
    :type manifold: bool
    """

    # Reload the injector from disk to ensure that this is starting with default values
    # This does only return an injector, but does not change the injector stored in config.injector
    newInjector = config.reload_injector()

    # Change the injector geometry. This is done for 'Oxygen' and 'Fuel' separately because the objects name differ,
    # however the actual change in the elements is the same.
    # The double asterisk ** is a notation for pow(), therefore change ** 2 is equal to pow(change, 2)
    if fluid == 'Oxygen':
        for i, element in enumerate(newInjector.LOX.elements):

            if manifold and element.type == 'Manifold':
                newInjector.LOX.elements[i].volume *= change

            if line_length and element.type == 'Line':
                newInjector.LOX.elements[i].length *= change

            if line_diameter and element.type == 'Line':
                newInjector.LOX.elements[i].diameter *= change
                newInjector.LOX.elements[i].area *= change ** 2

    if fluid == 'Fuel':
        for i, element in enumerate(newInjector.Fuel.elements):

            if manifold and element.type == 'Manifold':
                newInjector.Fuel.elements[i].volume *= change

            if line_length and element.type == 'Line':
                newInjector.Fuel.elements[i].length *= change

            if line_diameter and element.type == 'Line':
                newInjector.Fuel.elements[i].diameter += change
                newInjector.Fuel.elements[i].area += change ** 2

    # Recalculate all injector parameters to account for the changed geometry
    newInjector.recalculate_injector()

    # Compared to the previous code, there's a logical change regarding config.casMode. In the previous code, this
    # this parameter was set to True, later to False, then to True again (etc.). This is here not done at all.
    # However, the actual evaluation of config.casMode should return the same result.

    return newInjector

