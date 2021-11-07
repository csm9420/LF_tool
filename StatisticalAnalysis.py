# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 13:04:39 2020

@author: micha
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import pathlib
import csv


def write_unstable_points_to_file(path, unstable_points, N, tau_O, tau_F, tau_s):
    """ Writes all points in unstablePoints in the corresponding result file.

    """

    # New version - this does something different than the old code
    # Instead of writing multiple files, this creates just one single file containing all data
    # Multiply the time lag values by 1000 to convert to milliseconds

    path = os.path.join(path, 'Results.csv')
    header = ['N', 'tau_s', 'tau_Ox', 'tau_Fu', 'Frequency', 'Gain']
    if not os.path.exists(path):
        with open(path, 'a+', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=header, delimiter='\t')
            writer.writeheader()

    with open(path, 'a+', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=header, delimiter='\t')
        for frequency, gain in unstable_points:
            writer.writerow({'N': '{:.2f}'.format(N),
                             'tau_s': '{:.2f}'.format(tau_s * 1000),
                             'tau_Ox': '{:.2f}'.format(tau_O * 1000),
                             'tau_Fu': '{:.2f}'.format(tau_F * 1000),
                             'Frequency': '{:.2f}'.format(frequency),
                             'Gain': '{:.2f}'.format(gain)})

    #
    #
    # Creates the needed path (directories ans subdirectories)
    # path = os.path.join(path, 'FrequencyFiles')
    # filename = 'InteractionIndex={}_ResultFile.dat'.format(str(np.round(N, 2)))
    #
    # # If does not exists, create and write the header
    # if not os.path.exists(path):
    #     pathlib.Path(path).mkdir(parents=True, exist_ok=True)
    #     with open(os.path.join(path, filename), 'a+') as file:
    #         file.write("tauO\ttauF\ttauS\tN\tFrequency\tGain\n")
    #
    # # Write the data
    # with open(os.path.join(path, filename), 'a+') as file:
    #     for (frequency, gain) in unstable_points:
    #         file.write(
    #             str(format(np.round(tau_O, 5), ".5f")) + "\t" +
    #             str(format(np.round(tau_F, 5), ".5f")) + "\t" +
    #             str(format(np.round(tau_s, 5), ".5f")) + "\t" +
    #             str(format(np.round(N, 2), ".2f")) + "\t" +
    #             str(format(np.round(frequency, 2), ".2f")) + "\t" +
    #             str(format(np.round(gain, 2), ".2f")) + "\n")

    #
    #
    # workDir = os.getcwd()
    # path = workDir + "\\StabilityAnalysis"
    # if not (os.path.exists(path)):
    #     os.mkdir(path)
    # path = path + "\\"  + config.identifier
    # if not (os.path.exists(path)):
    #     os.mkdir(path)
    # path = path + "\\FrequencyFiles"
    # if not (os.path.exists(path)):
    #     os.mkdir(path)
    #
    # frequency_file = 'InteractionIndex={}_ResultFile.dat'.format(str(np.round(N, 2)))
    #
    # freqTxtPath = path + "\\" +"InteractionIndex=" + str(format(np.round(N, 2), '.2f')) + "_ResultFile.dat"
    #
    # Opening and writing header for the resultFile
    # if(os.path.exists(freqTxtPath)):
    #     freqTxt = open(freqTxtPath, "a")
    # else:
    #     try:
    #         freqTxt = open(freqTxtPath, "wt")
    #         freqTxt.write("tauO\ttauF\ttauS\tInteraction Index\tFrequency\tGain\n")
    #     except:
    #         freqTxt = open(freqTxtPath, "w")
    #
    # for point in unstablePoints:
    #     freqTxt.write(str(format(np.round(tau_O, 5), ".5f")) + "\t" + str(format(np.round(tau_F, 5), ".5f")) + "\t" + str(format(np.round(tau_s, 5), ".5f")) \
    #                   + "\t" + str(format(np.round(N, 2), ".2f")) + "\t" + str(format(np.round(point[0], 2), ".2f")) + "\t" + str(format(np.round(point[1], 2), ".2f")) + "\n")
    # freqTxt.close()


def createStatisticalPlot(path, freq_range, freq_tests, unstableArray, resolution):
    relativeAccumulation = []

    # for i in range(0, int(np.round((config.freq_range[1]) / resolution, 1))):
    for i in range(0, int(np.round((freq_range[1]) / resolution, 1))):
        relativeAccumulation.append(0)

    for i in range(len(unstableArray)):
        for j in range(len(unstableArray[i])):
            for k in range(len(unstableArray[i][j])):
                for l in range(len(unstableArray[i][j][k])):
                    for m in unstableArray[i][j][k][l]:
                        try:
                            slot = int(m[0] / resolution)
                            relativeAccumulation[slot] += m[1]
                        except:
                            for point in m:
                                slot = int(point[0] / resolution)
                                relativeAccumulation[slot] += point[1]

    fig = plt.figure()
    # plt.bar(np.arange(0, config.freq_range[1], resolution), relativeAccumulation, width=resolution * 0.9, align='center')
    plt.bar(np.arange(0, freq_range[1], resolution), relativeAccumulation, width=resolution * 0.9,
            align='center')
    plt.xlim(freq_range[0], freq_range[1])
    plt.ylim(0, max(relativeAccumulation) * 1.2)
    plt.grid()
    for m in freq_tests:
        plt.axvline(x=m, color='r')
    # plt.savefig("StabilityAnalysis\\" + config.identifier + "\\HistogramOverall" + str(resolution) + "Hz.png")
    plt.savefig(os.path.join(path, 'HistogramOverall{}Hz.png'.format(resolution)))
    plt.close(fig)


# Plot prints if a combination of tau_o and tau_f hit the test frequency
def hit_point_plot(path, hit_matrix, N):
    """ Create a plot tau_O over tau_F showing if a combination (tau_F, tau_O) hit the test frequency

    - Black point: Unstable in simulation
    - Red square: Unstable in simulation and tests
    - Green triangle: Stable

    :param path: Main result path
    :param hit_matrix:
    :param N: Interaction index
    :return:
    """
    x0 = []
    y0 = []
    x1 = []
    y1 = []
    x2 = []
    y2 = []

    tauFVal = [i[0][0] * 1000 for i in hit_matrix]
    tauOVal = [j[1] * 1000 for j in hit_matrix[0]]

    for i in hit_matrix:
        for j in i:
            if j[2] == 1:
                x1.append(j[0] * 1000)
                y1.append(j[1] * 1000)
            elif j[2] == 2:
                x2.append(j[0] * 1000)
                y2.append(j[1] * 1000)
            else:
                x0.append(j[0] * 1000)
                y0.append(j[1] * 1000)

    fig = plt.figure()
    # plt.grid(zorder=1)
    l1 = plt.scatter(x0, y0, s=25, c='k', zorder=2)
    l2 = plt.scatter(x1, y1, s=25, c='r', marker="s", zorder=2)
    l3 = plt.scatter(x2, y2, s=25, c='g', marker="^", zorder=2)
    plt.xlabel("tau_Fu\u2009[ms]")
    plt.ylabel("tau_Ox\u2009[ms]")
    plt.grid()
    # plt.xticks(tauFVal)
    # plt.yticks(tauOVal)
    leg = plt.legend([l1, l2, l3], ['unstable in simulation', 'unstable in simulation and test', 'stable'],
                     loc='lower center', ncol=3, bbox_to_anchor=(0.5, 1.01)).get_frame().set_linewidth(0.0)

    # Create the path, if it does not exists, but don't raise an error if the path already exists
    path = os.path.join(path, 'HitPlots')
    filename = 'HitPlot_N={}.png'.format(N)
    pathlib.Path(path).mkdir(parents=True, exist_ok=True)

    plt.savefig(os.path.join(path, filename))
    plt.close(fig)

    # workDir = os.getcwd()
    # path = workDir + "\\StabilityAnalysis"
    # if not (os.path.exists(path)):
    #     os.mkdir(path)
    # path = workDir + "\\StabilityAnalysis\\" + config.identifier
    # if not (os.path.exists(path)):
    #     os.mkdir(path)
    # path = path + "\\HitPlots"
    # if not (os.path.exists(path)):
    #     os.mkdir(path)
    # plt.savefig(path + "\\HitPlot" + str(N) + ".png")
    # plt.close(fig)
