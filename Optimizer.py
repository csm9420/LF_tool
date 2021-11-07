# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 12:37:00 2020

@author: micha
"""

import numpy as np
import NyquistKrit as nyq

from scipy.optimize import minimize


# Function optimizes the parameters tau_o, tau_f, tau_s and n for the minimal distance to known instability frequencies
# fRange:       vector of investigated frequency range
# Eq:           lambdified characteristic equation
# bounds:       parameter boundaries (tau_o, tau_f, tau_s/tau_o, n)
# freqs:        unstable frequencies found in experiment
def optimizeForFreq(fRange, Eq, boundaries, freqs):
    omegaVec = 2 * np.pi * fRange
    sVec = 1j * omegaVec

    def func(x):
        solVec = Eq(sVec, x[0], x[1], x[2], x[3], omegaVec)
        unstablePoints = nyq.open_loop_stability_criterion(solVec, fRange)
        dist = fRange[-1] - fRange[0]  # max distance possible
        for i in unstablePoints:
            for j in freqs:
                if (abs(i[0] - j) < dist):
                    dist = abs(i[0] - j)

        return dist

    startVal = [0.001, 0.001, 0.001, 1]
    cons = {'type': 'ineq', 'fun': func}

    res = minimize(func, startVal, bounds=boundaries)  # , constraints = cons)

    if (res.success):
        print(res)
    else:
        raise ValueError(res.message)

    print(func(res.x))

