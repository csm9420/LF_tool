# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 13:05:00 2020

@author: geuk_mi
"""

import os, sys


def ceaIdealCombustionTemp(fuel, T_iF, h_F, T_iO, h_O, rof, p_cc, accat):
    # Writing the input file
    os.chdir("CEAeinbindung")

    inputFile = open("cea.inp", "w")
    inputFile.write("problem\n")
    inputFile.write("case=NBK3 ro eq fac\n")
    inputString = "p(bar)=" + str(p_cc) + " o/f=" + str(rof) + " ac=" + str(accat) + "\n"
    inputFile.write(inputString)
    inputString = "react ox=O2 (L)   moles=1. O 2  t(k)=" + str(T_iO) + " h(j)=" + str(h_O) + "\n"
    inputFile.write(inputString)

    if (fuel == "CH4"):
        inputString = "      fu=" + fuel + " moles=1. C 1 H 4  t(k)=" + str(T_iF) + " h(j)=" + str(h_F) + "\n"
        inputFile.write(inputString)
        inputFile.write("only CH3 CH4 CH3OH CO CO2 COOH H HCO H2 HCHO H2 HCHO,formaldehyd HCOOH H2O OH O2\n")
    elif (fuel == "H2"):
        inputString = "      fu=" + fuel + " moles=1. H 2  t(k)=" + str(T_iF) + " h(j)=" + str(h_F) + "\n"
        inputFile.write(inputString)
        inputFile.write("only H2 O2 H2O OH H O H2O2 HO2 O3\n")

    inputFile.write("output trace=1.e-10 plot t\n")
    inputFile.write("end\n")
    inputFile.close()

    # Executing CEA
    os.system("cea2.exe")

    resultFile = open("cea.plt", "r")
    idealTemp = resultFile.readline()
    idealTemp = float(idealTemp)

    os.chdir("..")

    return idealTemp


###################################################################################################################################

def ceaRealCharacteristics(fuel, T_iF, h_F, T_iO, h_O, T_cc, rof, p_cc, accat):
    # Writing the input file
    os.chdir("CEAeinbindung")

    inputFile = open("cea.inp", "w")
    inputFile.write("problem\n")
    inputFile.write("case=NBK3 ro eq\n")
    inputString = "p(bar)=" + str(p_cc) + " o/f=" + str(rof) + " ac=" + str(accat) + " t(k)=" + str(T_cc) + "\n"
    inputFile.write(inputString)
    inputString = "react ox=O2 (L)   moles=1. O 2  t(k)=" + str(T_iO) + " h(j)=" + str(h_O) + "\n"
    inputFile.write(inputString)

    if (fuel == "CH4"):
        inputString = "      fu=" + fuel + " moles=1. C 1 H 4  t(k)=" + str(T_iF) + " h(j)=" + str(h_F) + "\n"
        inputFile.write(inputString)
        inputFile.write("only CH3 CH4 CH3OH CO CO2 COOH H HCO H2 HCHO H2 HCHO,formaldehyd HCOOH H2O OH O2\n")
    elif (fuel == "H2"):
        inputString = "      fu=" + fuel + " moles=1. H 2  t(k)=" + str(T_iF) + " h(j)=" + str(h_F) + "\n"
        inputFile.write(inputString)
        inputFile.write("only H2 O2 H2O OH H O H2O2 HO2 O3\n")

    inputFile.write("output trace=1.e-10 plot rho, m, gam, son, mach\n")
    inputFile.write("end\n")
    inputFile.close()

    # Executing CEA
    os.system("cea2.exe")

    solutionVector = []
    resultFile = open("cea.plt", "r")
    line = resultFile.readline().split("  ")

    # Writing solutions to vector
    solutionVector.append(float(line[1]))
    solutionVector.append(float(line[2]))
    solutionVector.append(float(line[3]))
    solutionVector.append(float(line[4]))
    solutionVector.append(float(line[5]))

    os.chdir("..")

    return solutionVector


##########################################################################################################################################

def ceaIdealCombustionTempExtended(fuel, T_iF, h_F, T_iO, h_O, rof, p_cc, accat):
    # Writing the input file
    os.chdir("CEAeinbindung")

    inputFile = open("cea.inp", "w")
    inputFile.write("problem\n")
    inputFile.write("ro eq fac\n")
    inputString = "p(bar)=" + str(p_cc) + " o/f=" + str(rof) + " ac=" + str(accat) + "\n"
    inputFile.write(inputString)
    inputString = "react ox=O2 (L)   moles=1. O 2  t(k)=" + str(T_iO) + " h(j)=" + str(h_O) + "\n"
    inputFile.write(inputString)

    if (fuel == "CH4"):
        inputString = "      fu=" + fuel + " moles=1. C 1 H 4  t(k)=" + str(T_iF) + " h(j)=" + str(h_F) + "\n"
        inputFile.write(inputString)
        inputFile.write("only CH3 CH4 CH3OH CO CO2 COOH H HCO H2 HCHO H2 HCHO,formaldehyd HCOOH H2O OH O2\n")
    elif (fuel == "H2"):
        inputString = "      fu=" + fuel + " moles=1. H 2  t(k)=" + str(T_iF) + " h(j)=" + str(h_F) + "\n"
        inputFile.write(inputString)
        inputFile.write("only H2 O2 H2O OH H O H2O2 HO2 O3\n")

    inputFile.write("output trace=1.e-10 plot t\n")
    inputFile.write("end\n")
    inputFile.close()

    # Executing CEA
    os.system("cea2.exe")

    solutionVector = []
    resultFile = open("cea.plt", "r")

    # Writing solutions to vector
    solutionVector.append(float(resultFile.readline()))
    solutionVector.append(float(resultFile.readline()))
    solutionVector.append(float(resultFile.readline()))

    os.chdir("..")

    return solutionVector


################################################################################################################################################

def ceaGetCStar(fuel, T_iF, h_F, T_iO, h_O, rof, p_cc, accat):
    # Writing the input file
    os.chdir("CEAeinbindung")

    inputFile = open("cea.inp", "w")
    inputFile.write("problem\n")
    inputFile.write("ro eq fac\n")
    inputString = "p(bar)=" + str(p_cc) + " o/f=" + str(rof) + " ac=" + str(accat) + "\n"
    inputFile.write(inputString)
    inputString = "react ox=O2 (L)   moles=1. O 2  t(k)=" + str(T_iO) + " h(j)=" + str(h_O) + "\n"
    inputFile.write(inputString)

    if (fuel == "CH4"):
        inputString = "      fu=" + fuel + " moles=1. C 1 H 4  t(k)=" + str(T_iF) + " h(j)=" + str(h_F) + "\n"
        inputFile.write(inputString)
        inputFile.write("only CH3 CH4 CH3OH CO CO2 COOH H HCO H2 HCHO H2 HCHO,formaldehyd HCOOH H2O OH O2\n")
    elif (fuel == "H2"):
        inputString = "      fu=" + fuel + " moles=1. H 2  t(k)=" + str(T_iF) + " h(j)=" + str(h_F) + "\n"
        inputFile.write(inputString)
        inputFile.write("only H2 O2 H2O OH H O H2O2 HO2 O3\n")

    inputFile.write("output trace=1.e-10 plot t\n")
    inputFile.write("end\n")
    inputFile.close()

    # Executing CEA
    os.system("cea2.exe")

    resultFile = open("cea.out", "r")

    currentLine = resultFile.readline()
    while not ("CSTAR" in currentLine):
        currentLine = resultFile.readline()
    resultVec = currentLine.split(" ")
    for i in resultVec:
        try:
            result = float(i)
            break
        except:
            pass

    os.chdir("..")

    return result


################################################################################################################################################

def ceaGetChamberSpeedValues(fuel, T_iF, h_F, T_iO, h_O, rof, p_cc, accat):
    # Writing the input file
    os.chdir("CEAeinbindung")

    inputFile = open("cea.inp", "w")
    inputFile.write("problem\n")
    inputFile.write("ro eq fac\n")
    inputString = "p(bar)=" + str(p_cc) + " o/f=" + str(rof) + " ac=" + str(accat) + "\n"
    inputFile.write(inputString)
    inputString = "react ox=O2 (L)   moles=1. O 2  t(k)=" + str(T_iO) + " h(j)=" + str(h_O) + "\n"
    inputFile.write(inputString)

    if (fuel == "CH4"):
        inputString = "      fu=" + fuel + " moles=1. C 1 H 4  t(k)=" + str(T_iF) + " h(j)=" + str(h_F) + "\n"
        inputFile.write(inputString)
        inputFile.write("only CH3 CH4 CH3OH CO CO2 COOH H HCO H2 HCHO H2 HCHO,formaldehyd HCOOH H2O OH O2\n")
    elif (fuel == "H2"):
        inputString = "      fu=" + fuel + " moles=1. H 2  t(k)=" + str(T_iF) + " h(j)=" + str(h_F) + "\n"
        inputFile.write(inputString)
        inputFile.write("only H2 O2 H2O OH H O H2O2 HO2 O3\n")

    inputFile.write("output trace=1.e-10 plot t\n")
    inputFile.write("end\n")
    inputFile.close()

    # Executing CEA
    os.system("cea2.exe")

    resultFile = open("cea.out", "r")

    currentLine = resultFile.readline()
    while not ("SON VEL" in currentLine):
        currentLine = resultFile.readline()
    resultVec = currentLine.split(" ")
    result = []
    for i in resultVec:
        try:
            result.append(float(i))
        except:
            pass

    soundSpeed = (result[1] + result[
        2]) / 2  # Average Soundspeed as middle between of values at combustion end and throat

    while not ("MACH NUMBER" in currentLine):
        currentLine = resultFile.readline()
    resultVec = currentLine.split(" ")
    result = []
    for i in resultVec:
        try:
            result.append(float(i))
        except:
            pass

    machNumber = result[1]
    meanVelo = soundSpeed * machNumber

    os.chdir("..")

    return soundSpeed, machNumber, meanVelo


################################################################################################################################################

def ceaGetSpecificHeatRatio(fuel, T_iF, h_F, T_iO, h_O, rof, p_cc, accat):
    # Writing the input file
    os.chdir("CEAeinbindung")

    inputFile = open("cea.inp", "w")
    inputFile.write("problem\n")
    inputFile.write("ro eq fac\n")
    inputString = "p(bar)=" + str(p_cc) + " o/f=" + str(rof) + " ac=" + str(accat) + "\n"
    inputFile.write(inputString)
    inputString = "react ox=O2 (L)   moles=1. O 2  t(k)=" + str(T_iO) + " h(j)=" + str(h_O) + "\n"
    inputFile.write(inputString)

    if (fuel == "CH4"):
        inputString = "      fu=" + fuel + " moles=1. C 1 H 4  t(k)=" + str(T_iF) + " h(j)=" + str(h_F) + "\n"
        inputFile.write(inputString)
        inputFile.write("only CH3 CH4 CH3OH CO CO2 COOH H HCO H2 HCHO H2 HCHO,formaldehyd HCOOH H2O OH O2\n")
    elif (fuel == "H2"):
        inputString = "      fu=" + fuel + " moles=1. H 2  t(k)=" + str(T_iF) + " h(j)=" + str(h_F) + "\n"
        inputFile.write(inputString)
        inputFile.write("only H2 O2 H2O OH H O H2O2 HO2 O3\n")

    inputFile.write("output trace=1.e-10 plot t\n")
    inputFile.write("end\n")
    inputFile.close()

    # Executing CEA
    os.system("cea2.exe")

    resultFile = open("cea.out", "r")

    currentLine = resultFile.readline()
    while not ("GAMMAs" in currentLine):
        currentLine = resultFile.readline()
    resultVec = currentLine.split(" ")
    result = []
    for i in resultVec:
        try:
            result.append(float(i))
        except:
            pass

    gamma = (result[1] + result[
        2]) / 2  # Average specific heat ration between end of combustion and throat #Isentropic exponent

    os.chdir("..")

    return gamma


##########################################################################################################################################################

def ceaGetMoleMassProduct(fuel, T_iF, h_F, T_iO, h_O, rof, p_cc, accat, throat=False):
    # Writing the input file
    os.chdir("CEAeinbindung")

    inputFile = open("cea.inp", "w")
    inputFile.write("problem\n")
    inputFile.write("ro eq fac\n")
    inputString = "p(bar)=" + str(p_cc) + " o/f=" + str(rof) + " ac=" + str(accat) + "\n"
    inputFile.write(inputString)
    inputString = "react ox=O2 (L)   moles=1. O 2  t(k)=" + str(T_iO) + " h(j)=" + str(h_O) + "\n"
    inputFile.write(inputString)

    if (fuel == "CH4"):
        inputString = "      fu=" + fuel + " moles=1. C 1 H 4  t(k)=" + str(T_iF) + " h(j)=" + str(h_F) + "\n"
        inputFile.write(inputString)
        inputFile.write("only CH3 CH4 CH3OH CO CO2 COOH H HCO H2 HCHO H2 HCHO,formaldehyd HCOOH H2O OH O2\n")
    elif (fuel == "H2"):
        inputString = "      fu=" + fuel + " moles=1. H 2  t(k)=" + str(T_iF) + " h(j)=" + str(h_F) + "\n"
        inputFile.write(inputString)
        inputFile.write("only H2 O2 H2O OH H O H2O2 HO2 O3\n")

    inputFile.write("output trace=1.e-10 plot t\n")
    inputFile.write("end\n")
    inputFile.close()

    # Executing CEA
    os.system("cea2.exe")

    resultFile = open("cea.out", "r")

    currentLine = resultFile.readline()
    while not ("M," in currentLine):
        currentLine = resultFile.readline()
    resultVec = currentLine.split(" ")
    for i in resultVec:
        try:
            result = float(i)
            if not (throat):
                break
        except:
            pass

    os.chdir("..")

    return result