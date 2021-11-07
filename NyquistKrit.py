# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 13:07:58 2020

@author: geuk_mi
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import Functions as func
from scipy.signal import find_peaks
import pathlib
from numpy import savetxt

from sympy import symbols, lambdify, Function, exp


# def NyquistDiagramm(config, G, tau_o, tau_f):
#     Re = [np.real(i) if not np.isnan(np.real(i)) else 0 for i in G]
#     Im = [np.imag(i) if not np.isnan(np.imag(i)) else 0 for i in G]
#     scale = max((max(np.abs(Re)), max(np.abs(Im))))
#     if (scale < 1.5):
#         scale = 1.5
#
#     fig = plt.figure()
#     fig.set_size_inches(10, 10)
#     plt.xlim(-scale, scale)
#     plt.ylim(-scale, scale)
#     plt.plot(Re, Im)
#
#     # circular grid
#     alpha = np.arange(0, 2 * np.pi, 0.01)
#     unityCircle = np.cos(alpha) + 1j * np.sin(alpha)
#     Circle05 = 0.5 * (np.cos(alpha) + 1j * np.sin(alpha))
#     Circle15 = 1.5 * (np.cos(alpha) + 1j * np.sin(alpha))
#     Circle20 = 2.0 * (np.cos(alpha) + 1j * np.sin(alpha))
#     Circle25 = 2.5 * (np.cos(alpha) + 1j * np.sin(alpha))
#     Circle30 = 3.0 * (np.cos(alpha) + 1j * np.sin(alpha))
#     plt.plot(np.real(unityCircle), np.imag(unityCircle), c="k")
#     plt.plot(np.real(Circle05), np.imag(Circle05), '--', c="k", linewidth=1, alpha=0.6)
#     plt.plot(np.real(Circle15), np.imag(Circle15), '--', c="k", linewidth=1, alpha=0.6)
#     plt.plot(np.real(Circle20), np.imag(Circle20), '--', c="k", linewidth=1, alpha=0.6)
#     plt.plot(np.real(Circle25), np.imag(Circle25), '--', c="k", linewidth=1, alpha=0.6)
#     plt.plot(np.real(Circle30), np.imag(Circle30), '--', c="k", linewidth=1, alpha=0.6)
#     plt.axvline(0, linestyle='--', linewidth=1, alpha=0.6, c='k')
#     plt.axhline(0, linestyle='--', linewidth=1, alpha=0.6, c='k')
#
#     workDir = os.getcwd()
#     path = workDir + "\\StabilityAnalysis"
#     if not (os.path.exists(path)):
#         os.mkdir(path)
#     path = workDir + "\\StabilityAnalysis\\" + config.identifier
#     if not (os.path.exists(path)):
#         os.mkdir(path)
#     path = path + "\\NyquistDiagrams"
#     if not (os.path.exists(path)):
#         os.mkdir(path)
#     path = path + "\\NyquistDiagramTau_O_" + str(format(np.round(tau_o, 5), '.5f')) + "_Tau_F_" + str(
#         format(np.round(tau_f, 5), '.5f')) + ".png"
#
#     plt.savefig(path)
#     plt.close(fig)


def bode_diagram(path, freq_range, main_freq, hit_accuracy, ruemmler, full_transfer_function, G, freq, tau_o, tau_f,
                 tau_s, N,
                 numOfHits, unstablePoints=None, testFreqs=[]):
    """ Creates a Bode diagram.

    """
    phi = np.arctan2(np.imag(G), np.real(G))
    phi = phi * 180 / np.pi

    # plt.figure()

    # Align labels on y-axis
    fig, _ = plt.subplots(2, 1)
    fig.align_ylabels()

    sb1 = plt.subplot(2, 1, 1)
    sb1.plot(freq, abs(G))
    # sb1.set_xlim(config.freq_range[0], config.freq_range[1])
    sb1.set_xlim(min(freq_range), max(freq_range))
    yLim = [abs(i) for i in G if not np.isnan(i)]
    sb1.set_ylim(0, max(yLim) * 1.1)

    # if config.ruemmler and config.fullTransferFunc == False:
    if ruemmler is True and full_transfer_function is False:
        # Threshold of 0.5 as defined by Rümmler
        sb1.axhline(0.5, c='k', linestyle='--')
    else:
        sb1.axhline(1, c='k', linestyle='--')
    sb1.grid()
    # sb1.set_ylabel("Gain |G|")
    sb1.set_ylabel("Gain [-]")

    # plt.title(
    #     "Tau_Ox = " + str(format(np.round(tau_o, 5), '.5f')) + " , Tau_F = " + str(format(np.round(tau_f, 5), '.5f')) + \
    #     " , Tau_S = " + str(format(np.round(tau_s, 5), '.5f')) + " , InterIndex = " + str(
    #         format(np.round(N, 1), '.1f')))

    # The '\u2009' is the unicode character for a thin space, which creates a nice formatting of the plot title
    plt.title(
        'tau_Ox = {:2.2f}\u2009ms   tau_Fu = {:2.2f}\u2009ms   tau_s = {:2.2f}\u2009ms   N = {:1.2f}\u2009[-]'.format(
            tau_o * 1000, tau_f * 1000, tau_s * 1000, N))

    sb2 = plt.subplot(2, 1, 2)
    sb2.plot(freq, phi)
    sb2.set_xlim(min(freq_range), max(freq_range))
    sb2.grid()
    sb2.set_xlabel("Frequency [Hz]")
    sb2.set_ylabel("Phase [°]")

    # If there are unstable frequencies, mark each with red cross 'x' in the diagram
    if unstablePoints:
        for point in unstablePoints:
            sb1.plot(point[0], point[1], marker="x", c="r")

    # If there are test results available
    for tFreq in testFreqs:
        if tFreq == main_freq:
            sb1.axvline(x=tFreq, color='r')
            # sb1.axvspan(tFreq - 40, tFreq + 40, color='#2ca02c', alpha=0.5)
            sb1.axvspan(tFreq - hit_accuracy, tFreq + hit_accuracy, color='#2ca02c', alpha=0.5)

        else:
            sb1.axvline(x=tFreq, color='b')
            # sb1.axvspan(tFreq - 40, tFreq + 40, color='#99FFFF', alpha=0.5)
            sb1.axvspan(tFreq - hit_accuracy, tFreq + hit_accuracy, color='#99FFFF', alpha=0.5)

    # workDir = os.getcwd()
    # path = workDir + "\\StabilityAnalysis"
    # if not (os.path.exists(path)):
    #     os.mkdir(path)
    # path = workDir + "\\StabilityAnalysis\\" + config.identifier
    # if not (os.path.exists(path)):
    #     os.mkdir(path)
    # path = path + "\\BodeDiagrams"
    # if not (os.path.exists(path)):
    #     os.mkdir(path)
    #
    # if unstablePoints == None and numOfHits == 0:
    #     path = path + "\\Stable"
    #     if not (os.path.exists(path)):
    #         os.mkdir(path)
    # else:
    #     path = path + "\\" + str(numOfHits) + "-Hits"
    #     if not (os.path.exists(path)):
    #         os.mkdir(path)
    #
    # path = path + "\\InteractionIndex = " + str(format(np.round(InterIndex, 1), '.1f'))
    # if not (os.path.exists(path)):
    #     os.mkdir(path)
    # path = path + "\\Tau_Ox = " + str(format(np.round(tau_o, 5), '.5f'))
    # if not (os.path.exists(path)):
    #     os.mkdir(path)
    # path = path + "\\BodeDiagramTau_F_" + str(format(np.round(tau_f, 5), '.5f')) + "_Tau_S_" + str(format(np.round(tau_s, 5), '.5f')) + ".png"
    # plt.savefig(path)

    # Creating the path to store the Bode Diagram, but don't raise an error if the path already exists
    path = os.path.join(path, 'BodeDiagrams')
    if unstablePoints is None and numOfHits == 0:
        path = os.path.join(path, 'Stable')
    else:
        path = os.path.join(path, '{}-Hits'.format(numOfHits))
    path = os.path.join(path, 'InteractionIndex = {:.1f}'.format(N), 'Tau_Ox = {:.5f}'.format(tau_o))
    pathlib.Path(path).mkdir(parents=True, exist_ok=True)

    # Save the plot
    filename = 'BodeDiagram_Tau_F = {:.5f} Tau_s = {:.5f}.png'.format(tau_f, tau_s)
    plt.savefig(os.path.join(path, filename))

    # Save instability points in excel files for comparing effect of the length of the line element
    freq_sp = []
    G_sp = []
    phi_sp = []
    for i in range(2, len(freq)):
        abs_G = abs(G[i])
        if (phi[i] * phi[i - 1] < 0 and phi[i] > phi[i - 1]):
            freq_sp.append(freq[i])
            G_sp.append(abs_G)
            phi_sp.append(phi[i])

    freq_sp = np.array(freq_sp)
    G_sp = np.array(G_sp)
    phi_sp = np.array(phi_sp)
    instability_points = np.vstack((freq_sp, G_sp, phi_sp)).T
    filename3 = 'InstabilityPoints_Tau_F = {:.5f} Tau_s = {:.5f}.csv'.format(tau_f, tau_s)
    np.savetxt(os.path.join(path, filename3), instability_points, delimiter=',', fmt='%f')

    # Save the bode plot result in text file
    bode_plot_result = np.vstack((freq, abs(G), phi)).T
    filename2 = 'BodeDiagraminTables_Tau_F = {:.5f} Tau_s = {:.5f}.csv'.format(tau_f, tau_s)
    np.savetxt(os.path.join(path, filename2), bode_plot_result, delimiter=',', fmt='%f')

def open_loop_stability_criterion(G, freq):
    """ Open loop stability criterion to identify all unstable frequencies.

        An unstable frequency is found, if for an frequency f the gain in higher than unity at a phase of -180 deg.
        See master thesis Geuking, DLR-LA-RAK-HF-RP-100, chapter 2.3.2 and listing B.3

        :return: List of all unstable points, whereby point is (frequency, gain)
        :rtype: list
    """
    phase = np.arctan2(np.imag(G), np.real(G)) * 180 / np.pi
    pointsOfn180 = [i for i in range(1, len(phase)) if
                    not np.sign(phase[i]) == np.sign(phase[i - 1]) and phase[i] > 120]
    unstablePoints = [[freq[j], np.abs(G[j])] for j in pointsOfn180 if np.abs(G[j]) >= 1]

    # Print the result, just for debugging
    # for i in pointsOfn180:
    #     if np.abs(G[i]) >= 1:
    #         print('f = {:<8.2} g = {:<8.2} phase = {:<8.2} phase_-1 = {:<8.2}'
    #               .format(freq[i], np.abs(G[i]), float(phase[i]), float(phase[i-1])))

    return unstablePoints


def ruemmler_criterion(config, G, freq):
    """ Criterion defined by Rümmler to identify all unstable frequencies.

    """
    peaks, _ = find_peaks(np.abs(G))

    phase_plus = np.arctan2(np.imag(G[peaks + 1]), np.real(G[peaks + 1]))
    phase_minus = np.arctan2(np.imag(G[peaks - 1]), np.real(G[peaks - 1]))
    delta_phase = np.zeros(len(peaks))

    # Distinguishing between different delta cases to get the right phase jump
    for q in range(len(peaks)):
        if (np.sign(phase_plus[q]) == np.sign(phase_minus[q])):  # Both angles positive or both negative
            delta_phase[q] = phase_plus[q] - phase_minus[q]
        elif (np.sign(phase_plus[q]) > np.sign(phase_minus[q])):  # Phase_plus positive and Phase_minus negative
            delta_phase[q] = phase_plus[q] - phase_minus[q]
            if (delta_phase[q] > np.pi):
                delta_phase[q] = 2 * np.pi - delta_phase[q]
        elif (np.sign(phase_plus[q]) < np.sign(phase_minus[q])):  # Phase_plus negative and Phase_minus positive
            delta_phase[q] = phase_plus[q] - phase_minus[q]
            if (delta_phase[q] < -np.pi):
                delta_phase[q] = -2 * np.pi - delta_phase[q]

    pos = func.find(delta_phase, lambda x: x >= 0)

    if (config.fullTransferFunc):
        gain = func.find(np.abs(G[peaks[pos]]), lambda x: x >= 1)
    else:
        gain = func.find(np.abs(G[peaks[pos]]), lambda x: x >= 0.5)
    unstablePoints = []
    for k in gain:
        unstablePoints.append([freq[peaks[pos[k]]], np.abs(G[peaks[pos[k]]])])

    return unstablePoints


def count_unstable_frequency_hits(freq_tests, accuracy, unstable_points):
    """ Counts how many unstable test frequencies are found in the simulation.

    """
    hits = 0
    for (simulation_frequency, simulation_gain) in unstable_points:
        for test_frequency in freq_tests:
            if test_frequency - accuracy < simulation_frequency < test_frequency + accuracy:
                # Frequency hit found
                hits += 1
    return hits

    # hits = 0
    # for p in unstable_points:
    #     for testF in config.freq_tests:
    #         if p[0] > testF - config.hitAccuracy and p[0] < testF + config.hitAccuracy:
    #             hits += 1
    # return hits

# Not sure what this was/is? Probably old code, comment out to prevent it from doing something unexpected/bad
#
# if __name__ == "__main__":
#     freqVec = np.arange(-50, 50, 0.01)
#     freqVecPos = np.arange(0,100, 0.01)
#     sVec    = 1j * (2 * np.pi * freqVec)
#     sVec2    = 1j * (2 * np.pi * freqVecPos)
#
#     s = symbols("s")
#     P = Function("P")("s")
#     P =  10*(s+1)/(s+10)
#     Plam = lambdify(s, P, "numpy")
#
#     Res = Plam(sVec)
#     #NyquistDiagramm(Res)
#     Res2 = Plam(sVec2)
#     bode_diagram(Res2, freqVecPos)
