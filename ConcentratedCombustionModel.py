# -*- coding: utf-8 -*-
"""
Low Frequency Tool

"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from tqdm import tqdm
import CharacteristicEquation
import sympy as sym
import NyquistKrit as nyq
import StatisticalAnalysis as stat
from Config import ParseConfig
from point_analysis import main_frequency_hit
from temperature_variation import temperature_variation
from injector_geometry_variation import injector_geometry_variation
from pressure_drop_iteration import PressureDropIteration

# Filename of / path to the excel file containing all general information
FILE_CONFIG = 'Settings.xlsx'

# Main directory where all results will be stored; based on the current directory;
# As example: The current director is 'C:\Users', the results are stored in 'C:\Users\Results'
RESULT_DIRECTORY = 'Results'

# Settings
matplotlib.use("Agg")
matplotlib.rc('figure', max_open_warning=0)
np.seterr(invalid='ignore')
np.seterr(divide='ignore')


if __name__ == "__main__":

    print('Starting low frequency tool')

    # Close all previous plots, if there are any
    plt.close('all')

    # Load the excel file defined by FILE_CONFIG containing all general information;
    # then renaming the index (= row's); replacing a str with the corresponding True/False boolean value;
    # general_config is type pd.DataFrame; possible to access some data by general_config.loc[index][column];
    # then removing all rows where the index start with a '#' to get a pd.DataFrame without comment rows
    general_config = pd.read_excel(FILE_CONFIG, engine='openpyxl')
    general_config = general_config \
        .replace(['True', 'true'], True) \
        .replace(['False', 'false'], False) \
        .rename(index=lambda k: general_config.loc[k]['# Identifier'])
    general_config = general_config.drop(list(filter(lambda k: k.startswith('#'), general_config.index.tolist())),
                                         axis='rows')

    # TODO: Maybe change the list(float) to a np.array for better performance
    # Iterate over all cells (here: item) and find any array values, then convert them into type list(float)
    # An array in the excel document has to start with '[' and end with ']',
    # values can be separated by '; ' or ';' or ', ' or ',' or ' '; as example [3.14; 42.7] or [3.14,12.8,9.1]
    for index, column in [(index, column) for index in general_config.index for column in general_config]:
        item = general_config.loc[index, column]
        if type(item) is str and item.startswith('[') and item.endswith(']'):
            general_config.at[index, column] = list(np.float_(re.split(r'[;,] ?| ', item[1:-1])))

    # From general_config, collect all columns which should be processed, then call parse_config() for each;
    # then add them to a list configs_to_process
    # The name of columns which should not be processed should start with a '#', for example '# unit'
    configs_to_process = []
    for column in general_config:
        if not column.startswith('#'):
            print('Loading configuration for load point {}'.format(column))
            configs_to_process.append(ParseConfig(column, RESULT_DIRECTORY, FILE_CONFIG, general_config[column]))

    # Dropping all variables which are not needed any more. This is definitely not necessary and has no impact on
    # the result, but it makes debugging a lot easier because there are less variables in scope
    del general_config, column, index, item, FILE_CONFIG, RESULT_DIRECTORY

    # Inside this loop, the actual calculation is performed;
    # Here, config contains all data for a single load point (a single column in the configuration excel file)
    for config in configs_to_process:

        # Close all previous plots, if there are any
        plt.close('all')

        # Just for display purposes
        print('Calculating load point {}'.format(config.identifier))

        injector = config.injector
        chamber = config.combustion_chamber

        if config.stability_analysis_mode:

            # Defining the vectors to create evenly (linear) split values between the limits defined in the
            # configuration file. Adding the 0.1 * step_width to the upper limit ensured that np.arange includes the
            # upper limit. Alternatively, np.linspace() could have been used but this does not return a np.array, and
            # is therefore not suitable.
            # Dividing the time lag vectors (tau_o_vec, tau_f_vec) by 1000 converts to seconds;
            # here now: time lag in seconds [s], frequency in [Hz]
            # Variables, which are vectors are named with the subfix '_vec' for better readability
            tau_O_vec = np.arange(min(config.tau_Ox), max(config.tau_Ox) + 0.1 * config.dtau_Ox, config.dtau_Ox) / 1000
            tau_F_vec = np.arange(min(config.tau_Fu), max(config.tau_Fu) + 0.1 * config.dtau_Fu, config.dtau_Fu) / 1000
            # This is should be set to a const of 0.08 ms in the input file
            tau_s_vec = np.arange(min(config.tau_s), max(config.tau_s) + 0.1 * config.dtau_s, config.dtau_s) / 1000
            N_vec = np.arange(min(config.N), max(config.N) + config.dN * 0.1, config.dN)
            frequency_vec = np.arange(min(config.freq_range), max(config.freq_range) + 0.1 * config.df, config.df)
            omega_vec = 2 * np.pi * frequency_vec
            s_vec = 1j * omega_vec

            # Defining the equations and the symbolic variables
            # To ensure clarification, all symbolic variables are named with the subfix '_sym'
            eqs = CharacteristicEquation.Equations(injector, config, chamber)
            s_sym, tau_O_sym, tau_F_sym, tau_s_sym, N_sym, omega_sym = sym.symbols(
                "s_sym tau_O_sym tau_F_sym tau_s_sym N_sym omega_sym")
            symChar = eqs.getSystemAdmittance(s_sym, tau_O_sym, tau_F_sym, tau_s_sym, N_sym, omega_sym)
            charEq = sym.lambdify((s_sym, tau_O_sym, tau_F_sym, tau_s_sym, N_sym, omega_sym), symChar, "numpy")

            # Dropping the symbolic variables because they're no longer needed. This is not really necessary but helps
            # for clarification - the variables created by the for... loops are representing further a specific value,
            # not a symbolic variable!
            del s_sym, tau_O_sym, tau_F_sym, tau_s_sym, N_sym, omega_sym, symChar

            # Create an object to monitor the progress; this is just for display and monitoring purposes
            monitor_stability_analysis = tqdm(total=len(tau_O_vec) * len(tau_F_vec) * len(tau_s_vec) * len(N_vec),
                                              desc='Stability analysis for {}'.format(config.identifier), leave=True,
                                              bar_format='{desc:<65}| {bar:25}| {n_fmt}/{total_fmt}')

            # Stores all unstable points for a given parameter combination of N, tau_F, tau_O, tau_s
            # A 'point' refers to (unstable_frequency, gain_at_this_frequency)
            unstablePointMatrix = np.zeros((len(N_vec), len(tau_F_vec), len(tau_O_vec), len(tau_s_vec)),
                                           dtype=tuple)

            # Iterating over N, tau_F, tau_O and tau_s to calculate a solution for every combination of these parameters
            for iterate_N, N in enumerate(N_vec):

                # Collection all results necessary for the HitPlot
                # 0: Unstable, but no hit
                # 1: Unstable and main unstable frequency hit
                # 2: Stable
                hitPlotMatrix = np.zeros((len(tau_F_vec), len(tau_O_vec)), dtype=np.ndarray)

                for iterate_tau_F, tau_F in enumerate(tau_F_vec):
                    for iterate_tau_O, tau_O in enumerate(tau_O_vec):
                        for iterate_tau_s, tau_s in enumerate(tau_s_vec):

                            # Calculate the solution for the current combination of N, tau_F, tau_O and tau_s
                            solution_vec = charEq(s_vec, tau_O, tau_F, tau_s, N, omega_vec)

                            # Determine all unstable frequencies in the solution. If both, ruemmler_criterion and
                            # nyquist_criterion, are False, then unstable_points is not assigned and will crash; this
                            # is intentional to enforce a user input in the configuration. A 'point' refers to
                            # (unstable_frequency, gain_at_this_frequency)
                            if config.open_loop_stability:
                                unstable_points = nyq.open_loop_stability_criterion(solution_vec, frequency_vec)
                            elif config.ruemmler_stability_criterion:
                                unstable_points = nyq.ruemmler_criterion(config, solution_vec, frequency_vec)

                            # Calculate how many unstable test frequencies are found
                            # If no unstable test frequencies are specified, this will be zero
                            number_of_hits = nyq.count_unstable_frequency_hits(config.freq_tests, config.hitAccuracy,
                                                                               unstable_points)

                            # Writing all unstable points = (frequency, gain) to the Frequency File
                            stat.write_unstable_points_to_file(config.result_path, unstable_points, N, tau_O, tau_F,
                                                               tau_s)

                            unstablePointMatrix[iterate_N][iterate_tau_F][iterate_tau_O][iterate_tau_s] = unstable_points

                            # TODO: better logic to create the bode diagram

                            # Creating the bode diagram
                            nyq.bode_diagram(path=config.result_path,
                                             freq_range=config.freq_range, main_freq=config.mainFreq,
                                             ruemmler=config.ruemmler_stability_criterion,
                                             full_transfer_function=config.fullTransferFunc,
                                             G=solution_vec, freq=frequency_vec,
                                             tau_o=tau_O, tau_f=tau_F, tau_s=tau_s, N=N,
                                             numOfHits=number_of_hits,
                                             unstablePoints=unstable_points,
                                             hit_accuracy=config.hitAccuracy,
                                             testFreqs=config.freq_tests, )

                            if len(unstable_points) == 0:
                                # No unstable points found, stable
                                hitPlotMatrix[iterate_tau_F][iterate_tau_O] = [tau_F, tau_O, 2]
                            elif main_frequency_hit(config.mainFreq, config.hitAccuracy, unstable_points):
                                # Unstable, and main unstable frequency from tests is found
                                hitPlotMatrix[iterate_tau_F][iterate_tau_O] = [tau_F, tau_O, 1]
                            else:
                                # Unstable, but main unstable frequency from tests is not found
                                hitPlotMatrix[iterate_tau_F][iterate_tau_O] = [tau_F, tau_O, 0]


                            # No unstable frequencies found, therefore stable simulation and make bode diagram for
                            # stable points is True
                            # if config.showStablePoints and len(unstable_points) == 0:
                            #     nyq.bode_diagram(path=config.result_path,
                            #                      freq_range=config.freq_range, main_freq=config.mainFreq,
                            #                      ruemmler=config.ruemmler_stability_criterion,
                            #                      full_transfer_function=config.fullTransferFunc,
                            #                      G=solution_vec, freq=frequency_vec,
                            #                      tau_o=tau_O, tau_f=tau_F, tau_s=tau_s, N=N,
                            #                      numOfHits=number_of_hits,
                            #                      hit_accuracy=config.hitAccuracy, testFreqs=config.freq_tests)
                            #     hitPlotMatrix[iterate_tau_F][iterate_tau_O] = [tau_F, tau_O, 2]

                            # Main unstable frequency hit
                            # if main_frequency_hit(config.mainFreq, config.hitAccuracy, unstable_points):
                            #     nyq.bode_diagram(config.result_path,
                            #                      config.freq_range, config.mainFreq, config.hitAccuracy,
                            #                      config.ruemmler_stability_criterion,
                            #                      config.fullTransferFunc,
                            #                      solution_vec, frequency_vec,
                            #                      tau_O, tau_F, tau_s, N,
                            #                      number_of_hits,
                            #                      unstable_points,
                            #                      config.freq_tests)
                            #     hitPlotMatrix[iterate_tau_F][iterate_tau_O] = [tau_F, tau_O, 1]

                            # Unstable frequencies in simulation, but not in the tests
                            # if not len(unstable_points) == 0 and len(config.freq_tests) == 0:
                            #     nyq.bode_diagram(config.result_path,
                            #                      config.freq_range, config.mainFreq, config.hitAccuracy,
                            #                      config.ruemmler_stability_criterion,
                            #                      config.fullTransferFunc,
                            #                      solution_vec, frequency_vec,
                            #                      tau_O, tau_F, tau_s, N,
                            #                      number_of_hits,
                            #                      unstable_points,
                            #                      config.freq_tests)

                            # Update monitoring, just for display purposes
                            monitor_stability_analysis.update()

                # For each interaction index N, create a HitPlot figure
                stat.hit_point_plot(config.result_path, hitPlotMatrix, N)

            # Create the Histogram plot
            stat.createStatisticalPlot(config.result_path, config.freq_range, config.freq_tests, unstablePointMatrix, 20)

            # Finish monitoring, just for display purposes
            monitor_stability_analysis.close()

        # TODO: Works, but add input parameter to settings file
        if config.pressureDropIterationMode:
            # sRM.pressureDropIteration(config)
            PressureDropIteration(config, False)

        if config.tempVariationMode:
            # sRM.temperatureVariation(config)
            temperature_variation(config,
                                  range_factor_low=config.tv_range_factor_low,
                                  range_factor_high=config.tv_range_factor_high,
                                  steps=config.tv_steps)

        if config.geometryVariationMode:
            # TODO: Why reset? This seems scary! - not sure about this, result is fine
            #   Don't use the code below - this is an old library
            # Injector Geometry has to be resetted if in CasMode to get initial geometric parameters
            # if config.casMode:
            #     config.casMode = False
            #     # injector = injector_configuration(injector_file, config)
            #     injector = config.injector
            #     config.casMode = True
            # sRM.geometricVariation(config, injector)
            injector_geometry_variation(config,
                                        range_factor_low=config.igv_range_factor_low,
                                        range_factor_high=config.igv_range_factor_high,
                                        steps=config.igv_steps,
                                        vary_line_diameter=config.igv_vary_line_diameter,
                                        vary_line_length=config.igv_vary_line_length,
                                        vary_manifold=config.igv_vary_manifold)

    print('Finished.')