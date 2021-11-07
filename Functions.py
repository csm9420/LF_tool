# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 13:03:07 2020

@author: geuk_mi
"""

import math
import numpy as np


def plinedrop(length, diameter, rho, u_pipe, nu):
    Re = (u_pipe * diameter) / nu
    delta_p_L = 0

    # Old formulation
    if (Re <= 2320):
        delta_p_L = 0.5 * (6 / Re) * rho * u_pipe ** 2 * (length / diameter)
    elif ((Re > 2320) and (Re <= 10 ** 5)):
        delta_p_L = 0.5 * (0.3164 / Re ** 0.25) * rho * (u_pipe) ** 2 * (length / diameter)
    elif ((Re > 10 ** 5) and (Re <= 5 * 10 ** 7)):
        lambda_temp = (1 / (1.819 * math.log10(Re) - 1.64)) ** 2
        delta_p_L = 0.5 * lambda_temp * rho * (u_pipe) ** 2 * (length / diameter)
    else:
        print("Re out of range")

    # print("Old: " , delta_p_L)

    # Reformulation of line pressure drop with explicit friction coefficient equation after Papaevangelou et al
    e_rough = 0.05 * 10 ** (-3)
    lambda_explicit = (0.2479 - 0.0000947 * (7 - math.log10(Re)) ** 4) / (
        math.log10(e_rough / (3.615 * diameter) + 7.366 / Re ** 0.9142)) ** 2
    delta_p_L = 0.5 * lambda_explicit * rho * (u_pipe) ** 2 * (length / diameter)

    # print("New: ", delta_p_L)

    return delta_p_L


##############################################################################################

# calculates the pressure drop across the orifice
def pdrop(rho, u, zv):
    return 0.5 * rho * u ** 2 * zv


##############################################################################################

# calculates the zeta value for the pressure drops
# it compares the cross areas of different pipes
# equations from VDI-Wärmeatlas
# A1 > A2 -> contraction
# A2 < A3 -> Expansion

def zeta_p(delta, A1, A2):
    if (A1 < A2):
        if ((delta >= 0) and (delta < math.pi / 6)):
            phi = (delta / (math.pi / 2)) + math.sin(2 * delta)

        elif ((delta >= math.pi / 6) and (delta <= math.pi / 2)):
            phi = 1.25 - delta / (2 * math.pi)

        zeta_drop = phi * (1 - A1 / A2) ** 2
    elif (A1 > A2):
        phi = A2 / A1
        zeta_drop = 0.4 * (1 - phi) * (2 * delta / np.pi) ** (1.83 * (1 - phi) ** 0.4)

    return zeta_drop


###############################################################################################

def find(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]