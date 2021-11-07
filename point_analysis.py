# -*- coding: utf-8 -*-
"""
Created on Fri May 28 14:58:45 2021

@author: spah_mr
"""

import numpy as np
import CharacteristicEquation
import sympy as sym
import NyquistKrit as nyq


def main_frequency_hit(main_frequency, accuracy, unstable_points):
    """ Check if the main unstable frequency was hit.

    :param main_frequency: main unstable frequency
    :param accuracy:
    :param unstable_points: all unstable (frequency, gain) values found by the simulation
    :return: True, if the main frequency was hit within the defined accuracy, False otherwise
    :rtype: bool
    """
    for (frequency, gain) in unstable_points:
        if main_frequency - accuracy < frequency < main_frequency + accuracy:
            return True

    # No hit, therefore return False
    return False


# Was prior present in front of the __main__, -> moved to extra file
# def run_point_analysis(config_it, injector_it=None):
def run_point_analysis(config, injector, chamber):
    # if config_it.geometryVariationMode:
    #     injector = injector_it
    # else:
    #     loc = ("Geometry_Injection_Configurations" + "\\" + config_it.injector_configuration + ".xlsx")
    #     wb = xlrd.open_workbook(loc)
    #     injector_file = wb.sheet_by_index(0)
    #     injector = injector_configuration(injector_file, config_it)
    #
    # chamber = CombustionChamber(config_it)
    #
    # df = config.df
    # fmin = config.freq_range[0]
    # fmax = config.freq_range[1]
    # freqs = np.arange(fmin, fmax + df * 0.1, df)
    # omegaVec = 2 * np.pi * freqs
    # sVec = 1j * omegaVec
    #
    # eqs = CharacteristicEquation.Equations(injector, config, chamber)
    # s, tau_o, tau_f, tau_s, n, omega = sym.symbols("s tau_o tau_f tau_s n omega")
    # symChar = eqs.getSystemAdmittance(s, tau_o, tau_f, tau_s, n, omega)
    # charEq = sym.lambdify((s, tau_o, tau_f, tau_s, n, omega), symChar, "numpy")
    #
    # tau_o_Vec = np.arange(config.tau_Ox[0], config.tau_Ox[1] + 0.1 * config.dtau_Ox, config.dtau_Ox) / 1000
    # tau_f_Vec = np.arange(config.tau_Fu[0], config.tau_Fu[1] + 0.1 * config.dtau_Fu, config.dtau_Fu) / 1000
    # tau_s_factors = np.arange(0.08, 0.081, 1) / 1000
    # nVec = np.arange(config.N[0], config.N[1] + config.dN * 0.01, config.dN)
    #
    # unstablePointMatrix = np.zeros((len(nVec), len(tau_f_Vec), len(tau_o_Vec), len(tau_s_factors)), dtype=tuple)
    # yCounter = 0
    # xCounter = 0
    # zCounter = 0
    # nCounter = 0

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

    # Stores how many of the unstable frequencies found in the tests are predicted by the simulation for a given
    # combination of N, tau_F, tau_O, tau_s
    # hitMatrix = np.zeros((len(N_vec), len(tau_F_vec), len(tau_O_vec), len(tau_s_value)))

    unstablePointMatrix = np.zeros((len(N_vec), len(tau_F_vec), len(tau_O_vec), len(tau_s_vec)),
                                   dtype=tuple)

    # for n in nVec:
    #     for tauF in tau_f_Vec:
    #         for tauO in tau_o_Vec:
    #             for tauS in tau_s_factors:
    #
    #                 solVec = np.zeros((len(sVec)), dtype=complex)
    #                 solVec = charEq(sVec, tauO, tauF, tauS, n, omegaVec)
    #
    #                 if config.ruemmler_stability_criterion:
    #                     unstablePoints = nyq.ruemmler_criterion(solVec, freqs)
    #                 else:
    #                     unstablePoints = nyq.open_loop_stability_criterion(solVec, freqs)
    #
    #                 # print(len(unstablePoints))
    #                 unstablePointMatrix[nCounter][xCounter][yCounter][zCounter] = unstablePoints
    #
    #                 zCounter += 1
    #             yCounter += 1
    #             zCounter = 0
    #         xCounter += 1
    #         yCounter = 0
    #         zCounter = 0
    #     nCounter += 1
    #     xCounter = 0
    #     yCounter = 0
    #     zCounter = 0

    # Iterating over N, tau_F, tau_O and tau_s to calculate a solution for every combination of these parameters
    for iterate_N, N in enumerate(N_vec):
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

                    unstablePointMatrix[iterate_N][iterate_tau_F][iterate_tau_O][iterate_tau_s] = unstable_points

    # if (config.geometryVariationMode):
    #     return unstablePointMatrix, solVec
    #
    # return unstablePointMatrix

    return unstablePointMatrix, solution_vec
