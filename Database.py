# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 13:12:00 2020

@author: geuk_mi
"""

import sqlite3
from sqlite3 import Error
import sys, os
import xml.etree.ElementTree as ET

currentTableName = ""


def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)
    return conn


#############################################################################################################################################################################

def executeSQLQuery(database, query):
    try:
        cursor = database.cursor()
        cursor.execute(query)
    except Error as e:
        print(e)
    return cursor.fetchall()


#############################################################################################################################################################################

def checkWritingTautToDatabase(database, taut_val):
    # overwriteMode = False
    def overwriteCheck():
        print("The Tau_t value: ", taut_val,
              " is already part of the database, do you want to overwrite or delete current data or ignore new data? (O/D/i)")
        inputVal = input()
        if (inputVal == 'O'):
            # overwriteMode = true
            return True
        elif (inputVal == 'D'):
            executeSQLQuery(database, "DELETE from " + currentTableName + " WHERE Taut = " + str(taut_val) + ";")
            return True
        elif (inputVal == 'i'):
            print("Current data will not be saved to database")
            return False
        else:
            print("Invalid input.")
            overwriteCheck()

    query = "SELECT * from " + currentTableName + " WHERE Taut = " + str(taut_val) + ";"
    if not (len(executeSQLQuery(database, query)) == 0):
        if (overwriteCheck()):
            return True
        else:
            return False

    return True


#############################################################################################################################################################################

def readPointsFromDataBase(database, conditionType, conditionValue, returnType1, returnType2=None, returnType3=None,
                           returnType4=None, valRange=False, conditionValEnd=None):
    if (conditionType == 'Tau_t'):
        conditionType = 'Taut'
    elif (conditionType == 'Tau_s'):
        conditionType = 'Taus'
    elif (conditionType == 'N'):
        conditionType = 'InteractionIndex'
    elif (conditionType == 'Freq'):
        conditionType = 'Frequency'
    elif (conditionType == 'Gain'):
        pass
    else:
        print("Invalid Condition Type. Exit program!")
        sys.exit()

    returnTypes = []

    if (returnType1 == 'Tau_t'):
        returnTypes.append('Taut')
    elif (returnType1 == 'Tau_s'):
        returnTypes.append('Taus')
    elif (returnType1 == 'N'):
        returnTypes.append('InteractionIndex')
    elif (returnType1 == 'Freq'):
        returnTypes.append('Frequency')
    elif (returnType1 == 'Gain'):
        returnTypes.append('Gain')
    else:
        print("Return type 1 is invalid. Exit program!")
        sys.exit()

    if not (returnType2 == None):
        if (returnType2 == 'Tau_t'):
            returnTypes.append('Taut')
        elif (returnType2 == 'Tau_s'):
            returnTypes.append('Taus')
        elif (returnType2 == 'N'):
            returnTypes.append('InteractionIndex')
        elif (returnType2 == 'Freq'):
            returnTypes.append('Frequency')
        elif (returnType2 == 'Gain'):
            returnTypes.append('Gain')
        else:
            print("Return type 2 is invalid. Exit program!")
            sys.exit()

    if not (returnType3 == None):
        if (returnType3 == 'Tau_t'):
            returnTypes.append('Taut')
        elif (returnType3 == 'Tau_s'):
            returnTypes.append('Taus')
        elif (returnType3 == 'N'):
            returnTypes.append('InteractionIndex')
        elif (returnType3 == 'Freq'):
            returnTypes.append('Frequency')
        elif (returnType3 == 'Gain'):
            returnTypes.append('Gain')
        else:
            print("Return type 3 is invalid. Exit program!")
            sys.exit()

    if not (returnType4 == None):
        if (returnType4 == 'Tau_t'):
            returnTypes.append('Taut')
        elif (returnType4 == 'Tau_s'):
            returnTypes.append('Taus')
        elif (returnType4 == 'N'):
            returnTypes.append('InteractionIndex')
        elif (returnType4 == 'Freq'):
            returnTypes.append('Frequency')
        elif (returnType4 == 'Gain'):
            returnTypes.append('Gain')
        else:
            print("Return type 4 is invalid. Exit program!")
            sys.exit()

    returnTypeQuery = ""
    for i in range(len(returnTypes)):
        returnTypeQuery += returnTypes[i]
        if (i < len(returnTypes) - 1):
            returnTypeQuery += ", "

    if not (valRange):
        query = "SELECT " + returnTypeQuery + ", " + conditionType + " from " + currentTableName + " WHERE " + conditionType + " = " + str(
            conditionValue) + ";"
        solValues = executeSQLQuery(database, query)
    else:
        query = "SELECT " + returnTypeQuery + ", " + conditionType + " from " + currentTableName + " WHERE " + conditionType + " BETWEEN " \
                + str(conditionValue) + " AND " + str(conditionValEnd) + ";"
        solValues = executeSQLQuery(database, query)

    resultVec = []
    for i in range(len(returnTypes) + 1):
        tmp = []
        for j in solValues:
            tmp.append(j[i])
        resultVec.append(tmp)

    return resultVec


#############################################################################################################################################################################

def customQuery(database, query, noOfReturns):
    solValues = executeSQLQuery(database, query)
    resultVec = []
    for i in range(noOfReturns):
        tmp = []
        for j in solValues:
            tmp.append(j[i])
        resultVec.append(tmp)

    return resultVec


#############################################################################################################################################################################

def writeDataSetToDatabase(database, tauO, tauF, tauS, n, freq, gain):
    if (checkWritingTautToDatabase(database, taut)):
        for i in range(len(tauS)):
            query = "INSERT INTO " + currentTableName + " (TauO, TauF, TauS, InteractionIndex, Frequency, Gain) VALUES (" \
                    + str(tauO) + "," + str(tauF) + "," + str(taus[i]) + "," + str(n[i]) + "," + str(
                freq[i]) + "," + str(gain[i]) + ");"
            executeSQLQuery(database, query)


#############################################################################################################################################################################

def getConfigFromXML(xmlName):
    tree = ET.parse(xmlName)
    root = tree.getroot()

    for elem in root:
        # Getting testCase name and removing iteration number at the end
        testCase = ""
        try:
            testCaseTmp = \
            [subsubelem.text for subsubelem in [subelem for subelem in elem if elem.tag == 'GeneralSettings'] if
             subsubelem.tag == 'RunName'][0]
        except:
            pass
        testCaseVec = testCaseTmp.split('_')
        testCaseVec.pop(-1)
        for i in range(len(testCaseVec)):
            testCase += testCaseVec[i]
            if i < len(testCaseVec) - 1:
                testCase += "_"

        # Getting Fuel name from xml
        try:
            fuelType = \
            [subsubelem.text for subsubelem in [subelem for subelem in elem if elem.tag == 'LoadPointDefinition'] if
             subsubelem.tag == 'Fuel'][0]
        except:
            pass

        # Getting injectorConfig name from xml
        try:
            injectorConfig = \
            [subsubelem.text for subsubelem in [subelem for subelem in elem if elem.tag == 'InjectorGeometry'] if
             subsubelem.tag == 'InjectorGeometryConfigFile'][0]
        except:
            pass

    if (fuelType == "CH4"):
        fuelName = "Methane"
    elif (fuelType == "H2"):
        fuelName = "Hydrogen"
    else:
        print("Invalid Fuel Type!")
        fuelName = ""

    return fuelName.replace('-', '_'), injectorConfig.replace('-', '_'), testCase.replace('-', '_')


#############################################################################################################################################################################

def openDatabaseAndTable(configFile):
    fuel, databaseName, tableName = getConfigFromXML(configFile)
    directory = os.path.dirname(__file__) + "\\StabilityAnalysis\\" + "BKD\\" + fuel
    db = create_connection(directory + "\\" + databaseName + ".db")
    global currentTableName
    currentTableName = tableName
    executeSQLQuery(db, """CREATE TABLE IF NOT EXISTS """ + tableName + """ ( Id integer PRIMARY KEY, Taut float NOT NULL, Taus float NOT NULL, InteractionIndex float NOT NULL, \
                                                                                    Frequency float NOT NULL, Gain float NOT NULL);""")
    return db


############################################################################################################################################################################

def getTableNames(database):
    query = "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;"

    return executeSQLQuery(database, query)


#############################################################################################################################################################################

def closeDatabase(database):
    database.commit()
    database.close()


#############################################################################################################################################################################

if __name__ == '__main__':
    db = openDatabaseAndTable("Configuration.xml")
    writeDataSetToDatabase(db, 0.5, 0.4, 1.2, 900, 2.0)
    closeDatabase(db)
