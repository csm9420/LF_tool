""" Stability research mode - temperature variation

Created on Fri Jun 3 14:03:21 2021

"""
import numpy as np
from point_analysis import run_point_analysis
import pathlib
import os
import matplotlib.pyplot as plt
from tqdm import tqdm


# Name of the directory where the results from the temperature variation will be stored
RESULT_DIRECTORY = 'TemperatureVariation'


def temperature_variation(config, range_factor_low, range_factor_high, steps):
    """ Calculates the stability analysis for a variation of the temperatures.

    This function varies the oxidizer injection temperature T_io and the fuel injection temperature between the default
    values given in config multiplied with the range factors.

    :param config:
    :type config: :class:`Config.ParseConfig`
    :param range_factor_low: Factor for the lower end of the temperature variation, e.g. 0.8
    :type range_factor_low: float
    :param range_factor_high: Factor for the higher limit of the temperature variation, e.g. 1.2
    :type range_factor_high: float
    :param steps: Amount of steps
    :type steps: int
    :rtype: None
    """

    new_config = config

    # Defining the temperature ranges for T_io and T_if, calculating the amount of steps by (high - low) / step
    # The conversion to int ensures that the input is valid, and values are rounded down => int(7.9) will evaluate to 7
    tempRangeOx = np.linspace(range_factor_low * new_config.T_io, range_factor_high * new_config.T_io, num=int(steps))
    tempRangeFu = np.linspace(range_factor_low * new_config.T_if, range_factor_high * new_config.T_if, num=int(steps))

    # Because the stability analysis should not be evaluated for multiple values of tau_Ox, tau_Fu and N, new values are
    # set to be the average of the input values
    new_config.tau_Ox = [(min(new_config.tau_Ox) + max(new_config.tau_Ox)) / 2, (min(new_config.tau_Ox) + max(new_config.tau_Ox)) / 2]
    new_config.dtau_Ox = 1
    new_config.tau_Fu = [(min(new_config.tau_Fu) + max(new_config.tau_Fu)) / 2, (min(new_config.tau_Fu) + max(new_config.tau_Fu)) / 2]
    new_config.dtau_Fu = 1
    new_config.N = [(min(new_config.N) + max(new_config.N)) / 2, (min(new_config.N) + max(new_config.N)) / 2]
    new_config.dN = 1

    # Define temperature range
    # startTempOx = 0.8 * new_config.T_io
    # endTempOx = 1.2 * new_config.T_io
    # tempRangeOx = np.linspace(startTempOx, endTempOx, 10)

    # startTempFu = 0.8 * new_config.T_if
    # endTempFu = 1.2 * new_config.T_if
    # tempRangeFu = np.linspace(startTempFu, endTempFu, 10)

    # Setting parameter values as average of input values
    # new_config.tau_Ox = [(new_config.tau_Ox[0] + new_config.tau_Ox[1]) / 2,
    #                          (new_config.tau_Ox[0] + new_config.tau_Ox[1]) / 2]
    # new_config.dtau_Ox = 1
    # new_config.tau_Fu = [(new_config.tau_Fu[0] + new_config.tau_Fu[1]) / 2,
    #                          (new_config.tau_Fu[0] + new_config.tau_Fu[1]) / 2]
    # new_config.dtau_Fu = 1
    # new_config.N = [(new_config.N[0] + new_config.N[1]) / 2, (new_config.N[0] + new_config.N[1]) / 2]
    # new_config.dN = 1

    # Create an object to monitor the progress. This is just for monitoring and display purposes
    monitor_progress = tqdm(total=steps * 2, desc='Temperature variation for {}'.format(new_config.identifier),
                            leave=True, bar_format='{desc:<65}| {bar:25}| {n_fmt}/{total_fmt}')

    # Calculate temperature variation
    pointsForTempOx = []

    for currentTempOx in tempRangeOx:
        # Calculate
        new_config.T_io = currentTempOx
        unstablePointsMat, _ = run_point_analysis(new_config, new_config.injector, new_config.combustion_chamber)

        # Assign results
        for i in range(len(unstablePointsMat)):
            for j in range(len(unstablePointsMat[i])):
                for k in range(len(unstablePointsMat[i][j])):
                    for l in range(len(unstablePointsMat[i][j][k])):
                        for m in unstablePointsMat[i][j][k][l]:
                            tmpPoint = [currentTempOx, m]
                            pointsForTempOx.append(tmpPoint)

        # Just for display purposes
        monitor_progress.update()

    pointsForTempFu = []
    for currentTempFu in tempRangeFu:
        # Calculate
        new_config.T_if = currentTempFu
        unstablePointsMat, _ = run_point_analysis(new_config, new_config.injector, new_config.combustion_chamber)

        # Assign results
        for i in range(len(unstablePointsMat)):
            for j in range(len(unstablePointsMat[i])):
                for k in range(len(unstablePointsMat[i][j])):
                    for l in range(len(unstablePointsMat[i][j][k])):
                        for m in unstablePointsMat[i][j][k][l]:
                            tmpPoint = [currentTempFu, m]
                            pointsForTempFu.append(tmpPoint)

        # Just for display purposes
        monitor_progress.update()

    # Finish monitoring, just for display purposes
    monitor_progress.close()

    # Create plots and save results

    # Create path to store results
    pathlib.Path(os.path.join(config.result_path, RESULT_DIRECTORY)).mkdir(parents=True)

    # Oxygen variation plot
    plt.figure()
    plt.grid()
    plt.xlim(0, 1000)
    plt.ylim(0, 5)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Gain Factor")
    x = [i[1][0] for i in pointsForTempOx]
    y = [i[1][1] for i in pointsForTempOx]
    col = [i[0] for i in pointsForTempOx]
    sc = plt.scatter(x, y, c=col, marker="+")
    cbar = plt.colorbar(sc)
    cbar.set_label("Oxygen Temperature [K]")
    plt.savefig(os.path.join(config.result_path, RESULT_DIRECTORY, 'OxygenTemperatureVariation'))

    # Fuel plot
    plt.figure()
    plt.grid()
    plt.xlim(0, 1000)
    plt.ylim(0, 5)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Gain Factor")
    x = [i[1][0] for i in pointsForTempFu]
    y = [i[1][1] for i in pointsForTempFu]
    col = [i[0] for i in pointsForTempFu]
    sc = plt.scatter(x, y, c=col, marker="+")
    cbar = plt.colorbar(sc)
    cbar.set_label("Fuel Temperature [K]")
    plt.savefig(os.path.join(config.result_path, RESULT_DIRECTORY, 'FuelTemperatureVariation'))
