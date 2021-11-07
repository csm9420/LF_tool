# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 14:30:00 2020

@author: micha
"""

import numpy as np
import sys

from sympy import sin, cos, exp, Function, sqrt, tanh


class Equations:

    def __init__(self, injector, config, chamber):

        # objects
        self.injector = injector
        self.config = config
        self.chamber = chamber

        # constants
        self.gc = 1  # 386.0874

        # Input for concentrated combustion
        self.Lc = self.chamber.length  # 60
        self.Dc = self.chamber.diameter  # 3.5
        self.Ac = np.pi / 4 * self.Dc ** 2
        self.pc = self.config.p_cc * 100000  # 1120
        self.mO = self.config.mDot_O2  # 6.64
        self.mF = self.config.mDot_F  # 9.06
        self.m = self.mO + self.mF

        # Additional parameters
        self.c = self.chamber.soundSpeed  # 5030*12
        self.u = self.chamber.meanVelocity  # 872*12
        self.rho = self.chamber.rho  # 0.00016725
        self.rhog = self.rho / self.gc
        self.gam = self.chamber.gamma  # 1.3778
        self.Me = self.u / self.c

        self.psi = 0.05 * self.Lc

        # ManifoldParameter
        self.rhoMO = self.injector.LOX.elements[0].rho  # 68.2 / 1728
        self.aMO = self.injector.LOX.elements[0].soundVelocity  # 2754 * 12
        self.VMO = self.injector.LOX.elements[0].volume  # 7.4
        self.KMO = self.aMO ** 2 * self.rhoMO / self.gc
        self.CmO = self.rhoMO * self.VMO / self.KMO
        #
        # TODO: Why is the velocity at the end of the chamber necessary for the manifold?
        self.MO = self.injector.LOX.elements[-1].velocity / self.aMO  # 52.1 * 12/self.aMO

        self.rhoMF = self.injector.Fuel.elements[0].rho  # 3.6 / 1728
        self.aMF = self.injector.Fuel.elements[0].soundVelocity  # 3559 * 12
        self.VMF = self.injector.Fuel.elements[0].volume  # 50.9
        self.KMF = self.aMF ** 2 * self.rhoMF / self.gc
        self.CmF = self.rhoMF * self.VMF / self.KMF
        self.MF = self.injector.Fuel.elements[-1].velocity / self.aMF  # 897 * 12/self.aMF

        # Injector Parameter
        self.pox = self.config.p_io * 100000  # 1264
        self.DpO = self.pox - self.pc
        self.AoO = self.injector.LOX.elements[-2].area  # self.injector.getAveragedArea("Oxygen")
        self.loO = self.injector.getLengthCas("Oxygen")

        self.RoO = 2 * self.DpO / self.mO
        self.LoO = self.loO / (self.AoO * self.gc)
        self.pf = self.config.p_if * 100000  # 1566
        self.DpF = self.pf - self.pc
        self.AoF = self.injector.Fuel.elements[-2].area  # self.injector.getAveragedArea("Fuel")
        self.loF = self.injector.getLengthCas("Fuel")
        self.RoF = 2 * self.DpF / self.mF
        self.LoF = self.loF / (self.AoF * self.gc)

        # Damping Parameters
    #        self.omega   = omega
    #        self.k       = self.omega / self.c
    #        self.lambda_ = 0 #2 * self.rho * 2554 * self.omega * self.omega ** (-1,19)

    def getSystemAdmittance(self, s, tau_o, tau_f, tau_s, n, omega):

        # Different modes for calculation of injector addmittance
        # CasMode: Single Manifold with "long orifice" after Dissertation by Casiano
        if (self.config.casMode):
            G_O = Function("G_O")(s)
            G_O = -1 / (self.RoO + self.LoO * s + 1 / (self.CmO * s))

            G_F = Function("G_F")(s)
            G_F = -1 / (self.RoF + self.LoF * s + 1 / (self.CmF * s))
        else:
            G_O = Function("G_O")(s)
            G_O = self.getInjectionSideTransferFunction(self.injector.LOX, s)

            G_F = Function("G_F")(s)
            G_F = self.getInjectionSideTransferFunction(self.injector.Fuel, s)

        # Careful!
        # These can be uncommented to examine only the oxygen injection side, or only the fuel injection side.
        #
        # G_O = 0
        # G_F = 0
        #

        self.Ge = ((self.gam + 1) / (2 * self.gam)) * (self.m / self.pc)
        self.Ae = self.Ac

        # Bi and Be are parameters for the modified Crocco model
        Bi = Function("Bi")(s)
        Bi = (((self.MO / self.aMO) * self.AoO * self.rhoMF + (
                    self.MF / self.aMF) * self.AoF * self.rhoMO - self.rhoMO * self.rhoMF / (
                           self.c * self.rho) * self.Ac) * self.gc \
              - self.rhoMF * G_O - self.rhoMO * G_F) \
             / (((self.MO / self.aMO) * self.AoO * self.rhoMF + (
                    self.MF / self.aMF) * self.AoF * self.rhoMO + self.rhoMO * self.rhoMF / (
                             self.c * self.rho) * self.Ac) * self.gc \
                - self.rhoMF * G_O - self.rhoMO * G_F)

        Be = Function("Be")(s)
        Be = (self.m * self.c / (self.gam * self.pc) - self.Ae * self.gc - self.Ge * self.c) / (
                    self.m * (self.c / (self.gam * self.pc)) + self.Ae * self.gc - self.Ge * self.c)

        GcMC = Function("GcMC")(s)
        GcMC = (((Be * exp(2 * s * (self.psi - self.Lc) / ((self.Me ** 2 - 1) * self.c)) + 1) / (
                    Be * exp(2 * s * (self.psi - self.Lc) / ((self.Me ** 2 - 1) * self.c)) - 1)) \
                - ((Bi * exp(2 * s * self.psi / ((self.Me ** 2 - 1) * self.c)) + 1) / (
                            Bi * exp(2 * s * self.psi / ((self.Me ** 2 - 1) * self.c)) - 1))) \
               * self.pc / (self.u * self.rhog * self.c) + self.pc / (self.rhog * self.c ** 2)

        # Hi and He are parameters for the damping model
        self.k = omega / self.c
        self.lambda_ = 2 * self.rho * 2554 * omega * omega ** (-1.19)

        self.aK = self.Me ** 2 - 1
        self.bK = - (3 / 2 * self.Me * self.lambda_ * self.k / (self.rho * s) + 2 * self.Me * self.k)
        self.cK = self.k ** 2 + self.k ** 2 * self.lambda_ / (s * self.rho)

        self.K1 = (-self.bK + sqrt(self.bK ** 2 - 4 * self.aK * self.cK)) / (2 * self.aK)
        self.K2 = (-self.bK - sqrt(self.bK ** 2 - 4 * self.aK * self.cK)) / (2 * self.aK)

        self.A2 = 1 / (self.rhog * self.c) * (self.lambda_ * self.Me / 2 + self.rho * self.K1 * s / self.k) / (
                    self.lambda_ + self.rho * (1 - self.Me * self.K1 / self.k) * s)
        self.B2 = 1 / (self.rhog * self.c) * (self.lambda_ * self.Me / 2 + self.rho * self.K2 * s / self.k) / (
                    self.lambda_ + self.rho * (1 - self.Me * self.K2 / self.k) * s)

        Hi = Function("Hi")(s)
        Hi = - ((((self.MO / self.aMO) * self.AoO * self.rhoMF + (
                    self.MF / self.aMF) * self.AoF * self.rhoMO + self.rhoMO * self.rhoMF / self.gc * self.Ac * self.B2) * self.gc \
                 - self.rhoMF * G_O - self.rhoMO * G_F) \
                / (((self.MO / self.aMO) * self.AoO * self.rhoMF + (
                            self.MF / self.aMF) * self.AoF * self.rhoMO + self.rhoMO * self.rhoMF / self.gc * self.Ac * self.A2) * self.gc \
                   - self.rhoMF * G_O - self.rhoMO * G_F))

        He = Function("He")(s)
        He = - (self.m / (self.gam * self.pc) + self.rho * self.Ae * self.B2 - self.Ge) / (
                    self.m / (self.gam * self.pc) + self.rho * self.Ae * self.A2 - self.Ge) \
             * exp(self.Lc / omega * s * (self.K1 - self.K2))

        GcDM = Function("GcDM")(s)
        GcDM = ((He * self.A2 * exp(- self.psi / omega * s * (self.K1 - self.K2)) + self.B2) / (
                    He * exp(- self.psi / omega * s * (self.K1 - self.K2)) + 1) \
                - (Hi * self.A2 * exp(- self.psi / omega * s * (self.K1 - self.K2)) + self.B2) / (
                            Hi * exp(- self.psi / omega * s * (self.K1 - self.K2)) + 1)) \
               * self.pc / self.u + self.pc / (self.rhog * self.c ** 2)

        Gj = Function("Gj")(s, tau_o, tau_f)
        Gj = (self.pc / self.m) * (G_O * exp(-s * tau_o) + G_F * exp(-s * tau_f))

        Gb = Function("Gb")(n, tau_s)
        Gb = n * (1 - exp(-s * tau_s))

        C_MC = Function("C_MC")
        C_MC = - (Gj + Gb) / GcMC

        C_DM = Function("C_DM")
        C_DM = - (Gj + Gb) / GcDM  # Open loop transfer function

        C_Ruem = Function("C_Ruem")
        if (self.config.fullTransferFunc):
            C_Ruem = - (Gj + Gb) / GcDM / (1 - ((Gj + Gb) / GcDM))  # Closed loop full transfer function
        else:
            C_Ruem = 1 / (1 - ((Gj + Gb) / GcDM))  # Closed loop only denominator

        if (self.config.ruemmler_stability_criterion):
            return C_Ruem
        else:
            return C_DM

    ##########################################################################################################

    def getInjectionSideTransferFunction(self, injectionSide, s):

        G = Function('G')(s)
        G = 0
        # New method using impedance approach (different from M. Geuking's method)
        G = (1) / self.getInjectionSideImpedance(injectionSide, s)
        return G

    ##########################################################################################################

    def getInjectionSideImpedance(self, injectionSide, s):
        Z_tot = Function('Z_tot')(s)
        if (injectionSide.type == "LOX"):
            massFlow = self.mO
        else:
            massFlow = self.mF

        if self.config.use_feed_system is True:
            Zt = Function('Zt')(s)
            Zt = 0
            Z_mani = Function('Z_mani')(s)
            for i in range(injectionSide.num_elements):
                element = injectionSide.elements[i]
                if (i < 6):
                    if (injectionSide.elements[i].type == "Line"):
                        Zt = self.getElementImpedance(element, s, Zt)
                    else:
                        Zt += self.getElementImpedance(element, s, Zt)
                elif (i < 7):
                    if (injectionSide.elements[6].type == "Manifold"):
                        # original
                        Zt = 1 / (1 / Zt + 1 / (self.getElementImpedance(element, s, Zt))) * injectionSide.elements[-1].numOfElements
                    else:
                        # Considering manifold as a low pass filter
                        Z_mani = 1 / ((1/(element.soundVelocity / (injectionSide.elements[i-2].area))) +
                                       (1/((-1) / ((element.length * element.area /
                                       pow(element.soundVelocity, 2)) * s)))) / 9.81
                        Zt = 9.81 * (Zt / 9.81 - Z_mani * tanh(s * element.length / element.soundVelocity)) \
                            / (1 - Zt / 9.81 / Z_mani * tanh(s * element.length / element.soundVelocity)) * injectionSide.elements[-1].numOfElements
                else:
                    if (injectionSide.elements[i].type == "Line"):
                        Zt = self.getElementImpedance(element, s, Zt)
                    else:
                        Zt += self.getElementImpedance(element, s, Zt)
            Z_tot = Zt/injectionSide.elements[-1].numOfElements
        else:
            Zt = Function('Zt')(s)
            Zt = 0
            switchCount = 0
            numParElements = 1
            for i in range(injectionSide.num_elements):
                element = injectionSide.elements[i]
                if (injectionSide.elements[i].type == "Line"):
                    Zt = self.getElementImpedance(element, s, Zt)
                else:
                    if (numParElements != injectionSide.elements[i].numOfElements):
                        Zt = Zt * injectionSide.elements[i].numOfElements
                        Zt += self.getElementImpedance(element, s, Zt)
                        numParElements = injectionSide.elements[i].numOfElements
                    else:
                        Zt += self.getElementImpedance(element, s, Zt)
            Z_tot = Zt / injectionSide.elements[-1].numOfElements

        return Z_tot

    ##########################################################################################################

    def getElementImpedance(self, element, s, Z_before):
        if (element.type == "Manifold"):
            Z = self.getManifoldImpedance(element, s)
        elif (element.type == "Orifice"):
            Z = self.getOrificeImpedance(element, s)
        elif (element.type == "Line"):
            Z = self.getLineImpedance(element, s, Z_before)
        else:
            return
        return Z

    def getManifoldImpedance(self, element, s):
        Z = (-1) / ((element.volume / pow(element.soundVelocity, 2)) * s)
        return Z

    def getOrificeImpedance(self, element, s):  # constant is (-2) in the reference of Casiano
        Z = (-1*2) * element.dp / element.mDot
        return Z

    def getLineImpedance(self, element, s, Z_before):
        Z0 = element.soundVelocity / element.area / 9.81
        Z_before = Z_before / 9.81
        Z = 9.81 * (Z_before - Z0 * tanh(s * element.length / element.soundVelocity))\
            /(1 - Z_before/Z0 * tanh(s * element.length / element.soundVelocity))
        return Z

    ##########################################################################################################