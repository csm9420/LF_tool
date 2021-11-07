# -*- coding: utf-8 -*-
"""
Created on Fri May 21 14:58:45 2021

@author: spah_mr
"""
import os
import sys
import yaml
import csv
import GUI
from Injector import injector_configuration, save_injector_elements
from Chamber import CombustionChamber
import numpy as np
import shutil
import pathlib
import time
import stat
from flow_solver import FlowSolverInjector

# Filename
CONFIG_FILE_NAME = 'Configurations'

# Directory to store the InputConfigurations
INPUT_CONFIGURATIONS_DIRECTORY = 'InputConfigurations'


class ParseConfig:
    """ Parsing all settings for a single load point.

    """

    def __init__(self, identifier, result_directory, settings_file, data):
        """ Creates a configuration.

        Creates a configuration element which contains all general settings, the injector and the combustion chamber
        for a single load point.

        :param identifier: Identifier, or name of the load point
        :type identifier: str
        :param data: Contains all information for this load point
        :type data: :class:pandas.Series
        :rtype : ParseConfig

        """

        self.identifier = identifier

        # Path to the folder where all results are stored
        self.result_path = os.path.join(os.getcwd(), result_directory, self.identifier)

        # Process
        self.stability_analysis_mode = data.loc['stability_analysis']
        self.pressureDropIterationMode = data.loc['pressureDropIterationMode']
        self.tempVariationMode = data.loc['tempVariationMode']
        self.geometryVariationMode = data.loc['geometryVariationMode']
        self.casMode = data.loc['casMode']

        # Overwrite old results
        self.overwrite_old_results = data.loc['overwrite']

        # Data evaluation
        self.fullTransferFunc = data.loc['fullTransferFunction']
        self.showStablePoints = data.loc['showStablePoints']
        self.hitAccuracy = data.loc['hitAccuracy']
        self.ruemmler_stability_criterion = data.loc['ruemmler_criterion']
        self.open_loop_stability = data.loc['open_loop_stability']

        # Previous values were hard-coded; now in excel-config file
        # self.ruemmler = False
        # self.fullTransferFunc = False
        # self.casMode = False
        # self.geometryVariationMode = False
        # self.tempVariationMode = False
        # self.showStablePoints = True
        # self.hitAccuracy = 10

        # General settings
        # self.recalcCharEq = data.loc['recalcCharEq']  # Recalculate characteristic equation
        self.freq_range = data.loc['freq_range']  # Frequency range [Hz]
        self.df = data.loc['df']  # Frequency Step [Hz]
        self.tau_s = data.loc['tau_s']  # Sensitive time lag [ms]
        self.dtau_s = data.loc['dtau_s']  # Sensitive time lag step [ms]
        self.tau_Ox = data.loc['tau_Ox']  # Time lag range oxidizer [ms]
        self.dtau_Ox = data.loc['dtau_Ox']  # Time lag step oxidizer [ms]
        self.tau_Fu = data.loc['tau_Fu']  # Time lag range fuel [ms]
        self.dtau_Fu = data.loc['dtau_Fu']  # Time lag step fuel [ms]
        self.N = data.loc['N']  # Interaction index [-]
        self.dN = data.loc['dN']  # Interaction index step [-]

        # Chamber geometry
        self.d_cc = data.loc['d_cc']  # Chamber diameter [mm]
        self.l_cc = data.loc['l_cc']  # Chamber length [mm]
        self.epsilon_cc = data.loc['epsilon_c']  # Contraction ratio [-]

        # Injector geometry
        self.injector_configuration = data.loc['injector_configuration']  # Injector geometry configuration file name

        # Load point definition
        self.fuel = data.loc['fuel']  # Fuel ('CH4' or 'H2')
        self.p_cc = data.loc['p_cc']  # Chamber pressure [bar]
        self.rof_cc = data.loc['rof_cc']  # ROF at primary pressure
        self.mDot_O2 = data.loc['mDot_O2']  # LOX mass flow [kg/s]
        self.mDot_F = data.loc['mDot_F']  # Fuel mass flow [kg/s]
        self.mDot_SF = data.loc['mDot_SF']  # Optional; secondary fuel mass flow [kg/s]
        self.fluid_wc = data.loc['fluid_wc']  # Window cooling fluid
        self.mDot_WC = data.loc['mDot_WC']  # Window cooling mass flow [kg/s]
        self.p_io = data.loc['p_io']  # Optional; LOX dome pressure [bar]
        self.p_if = data.loc['p_if']  # Optional; fuel dome pressure [bar]
        self.p_sf = data.loc['p_sf']  # ?
        self.p_wc = data.loc['p_wc']  # Optional; window cooling dome pressure [bar]
        self.T_io = data.loc['T_io']  # LOX injection temperature [K]
        self.T_if = data.loc['T_if']  # Fuel injection temperature [K]
        self.T_sf = data.loc['T_sf']  # ?
        self.T_wc = data.loc['T_wc']  # Window cooling injection temperature [K]
        self.eta_cc = data.loc['eta_cc']  # Combustion efficiency [-]

        # Frequency comparison
        self.mainFreq = data.loc['mainFreq']  # Main frequency found in tests [Hz]
        self.freq_tests = data.loc['freq_tests']  # Frequencies found in tests [Hz]

        # Temperature variation (tv)
        self.tv_range_factor_low = data.loc['tv_range_factor_low']  # Low range factor [-]
        self.tv_range_factor_high = data.loc['tv_range_factor_high']  # High range factor [-]
        self.tv_steps = data.loc['tv_steps']  # Amount of steps between low an high range factors [-]

        # Injector geometry variation (igv)
        self.igv_range_factor_low = data.loc['igv_range_factor_low']  # Low range factor [-]
        self.igv_range_factor_high = data.loc['igv_range_factor_high']  # High range factor [-]
        self.igv_steps = data.loc['igv_steps']  # Amount of steps between low and high range factors [-]
        self.igv_vary_line_length = data.loc['igv_vary_line_length']  # Determines if the line length should be varied
        self.igv_vary_line_diameter = data.loc[
            'igv_vary_line_diameter']  # Determines if the line diameter should be varied
        self.igv_vary_manifold = data.loc['igv_vary_manifold']  # Determines if the manifold volume should be varied

        # New, test parameters
        # self.p_ft = data.loc['p_ft']  # Pressure in the propellant/fuel tank [bar]
        # self.T_ft = data.loc['T_ft']  # Temperature in the propellant/fuel tank [K]
        # self.p_ot = data.loc['p_ot']  # Pressure in the oxidizer tank [bar]
        # self.T_ot = data.loc['T_ot']  # Temperature in the oxidizer tank [K]
        # self.injector_start_index = data.loc['injector_start_index']  # Index of the element which starts the injector
        # self.use_whole_feed_system = data.loc['use_whole_feed_system']

        self.use_flow_solver = data.loc['use_flow_solver']  # Determines if the flow solver should be used. If False, use Geuking's approach
        self.use_feed_system = data['use_feed_system']  # Determines whether the feed system should be used. If False, use only the injector part for the transfer function
        self.p_feed_o1 = data['p_feed_o1']  # Feed system oxygen manifold pressure [bar]
        self.p_feed_o2 = data['p_feed_o2']  # Feed system oxygen line connection pressure [bar]
        self.p_feed_f1 = data['p_feed_f1']  # Feed system fuel manifold pressure [bar]
        self.p_feed_f2 = data['p_feed_f2']  # Feed system fuel line connection pressure [bar]
        self.T_feed_o1 = data['T_feed_o1']  # Feed system oxygen manifold temperature [K]
        self.T_feed_o2 = data['T_feed_o2']  # Feed system oxygen line connection temperature [K]
        self.T_feed_f1 = data['T_feed_f1']  # Feed system fuel manifold temperature [K]
        self.T_feed_f2 = data['T_feed_f2']  # Feed system fuel line connection temperature [K]

        # TODO: Temporary fix because nyquist can't handle np.nan, expecting a list
        if np.isnan(np.sum(self.freq_tests)):
            self.freq_tests = []

        # To prevent the accidental overwriting of already existing solutions, check if there's any solution with the
        # same name already present; if yes - ask the user if it is okay to overwrite the data
        # If some needed paths do not exists, create them
        if os.path.exists(os.path.join(self.result_path)):
            if self.overwrite_old_results is not True:
                if not GUI.check_overwrite(self.identifier):
                    # The project already exists and overwriting data is not okay, therefore exit the whole software
                    sys.exit()
                else:
                    # Delete old files
                    shutil.rmtree(self.result_path, onerror=handle_file_deleting_error)
            else:
                # The path exists, and overwriting is okay, therefore, delete the whole path (directory including all
                # subdirectories and files)
                #
                # If there is any windows explorer open in the folder, rmtree will raise an error. This catches the
                # error and gives the Windows explorer enough time to reload itself automatically.
                try:
                    shutil.rmtree(self.result_path, onerror=handle_file_deleting_error)
                except WindowsError as e:
                    pass
                time.sleep(1)

                # Create the path to have a new, empty folder
                os.mkdir(self.result_path)
        else:
            # Path does not exists, therefore overwriting is not possible, but the directory(s) has to be created
            pathlib.Path(self.result_path).mkdir(parents=True)

        # Create the directory where the InputConfigurations should be stored
        pathlib.Path(os.path.join(self.result_path, INPUT_CONFIGURATIONS_DIRECTORY)).mkdir(parents=True)

        # Create the injector

        if self.use_flow_solver is True:
            path_flow_solver = os.path.join(self.result_path, 'flow_solver')
            pathlib.Path(path_flow_solver).mkdir(parents=True)
            self.injector = FlowSolverInjector(path=path_flow_solver, settings_file=settings_file,
                                               sheet=self.injector_configuration,
                                               mDot_LOX=self.mDot_O2, mDot_fuel=self.mDot_F,
                                               p_io=self.p_io * 1e5, p_if=self.p_if * 1e5, p_cc=self.p_cc * 1e5,
                                               p_feed_o1=self.p_feed_o1 * 1e5, p_feed_f1=self.p_feed_f1 * 1e5,
                                               p_feed_o2=self.p_feed_o2 * 1e5, p_feed_f2=self.p_feed_f2 * 1e5,
                                               T_feed_f1=self.T_feed_f1, T_feed_o1=self.T_feed_o1,
                                               T_feed_f2=self.T_feed_f2, T_feed_o2=self.T_feed_o2,
                                               T_io=self.T_io, T_if=self.T_if, fuel_type=self.fuel,
                                               p_sf=self.p_sf * 1e5, p_wc=self.p_wc * 1e5,
                                               T_sf=self.T_sf, T_wc=self.T_wc,
                                               fluid_wc=self.fluid_wc,
                                               mDot_sf=self.mDot_SF, mDot_wc=self.mDot_WC,
                                               use_feed_system=self.use_feed_system)
        else:
            # Writing all input parameters like this ensures that nothing is accidentally switched
            self.injector = injector_configuration(self, injector_filename=self.injector_configuration,
                                                   fluid_fuel=self.fuel,
                                                   fluid_window_cooling=self.fluid_wc,
                                                   T_io=self.T_io, T_if=self.T_if, T_sf=self.T_sf, T_wc=self.T_wc,
                                                   mDot_O2=self.mDot_O2, mDot_F=self.mDot_F, mDot_SF=self.mDot_SF,
                                                   mDot_WC=self.mDot_WC)
            # sys.exit()

        # Save all injector elements
        path_injector = os.path.join(self.result_path, 'Injector')
        pathlib.Path(path_injector).mkdir(parents=True)
        save_injector_elements(path_injector, 'LOX_Injector', self.injector.LOX.elements)
        save_injector_elements(path_injector, 'Fuel_Injector', self.injector.Fuel.elements)

        # sys.exit()

        # self.injector = self.reload_injector()

        # Create the combustion chamber
        # Writing all input parameters like this ensures that nothing is accidentally switched
        self.combustion_chamber = CombustionChamber(self, d_cc=self.d_cc, l_cc=self.l_cc, epsilon_cc=self.epsilon_cc,
                                                    p_cc=self.p_cc, eta_cc=self.eta_cc, rof=self.rof_cc,
                                                    mDot_O2=self.mDot_O2, mDot_fuel_primary=self.mDot_F, fuel=self.fuel,
                                                    T_if=self.T_if, T_io=self.T_io)

        # Write all data to a yaml file The str()-wrapping around the values is necessary
        # because PyYaml tries to use numpy-values for some attributes, but is unable to display these afterwards (
        # kinda stupid behavior of the PyYaml library!)
        with open(os.path.join(self.result_path, INPUT_CONFIGURATIONS_DIRECTORY, CONFIG_FILE_NAME + '.yaml'),
                  'w') as yaml_file:
            yaml_content = yaml.dump({
                'Identifier': str(self.identifier),
                'Process': {
                    'defaultMode': str(self.stability_analysis_mode),
                    'pressureDropIterationMode': str(self.pressureDropIterationMode),
                    'tempVariationMode': str(self.tempVariationMode),
                    'geometryVariationMode': str(self.geometryVariationMode),
                    'casMode': str(self.casMode)
                },
                'Data evaluation': {
                    'fullTransferFunction': str(self.fullTransferFunc),
                    'showStablePoints': str(self.showStablePoints),
                    'hitAccuracy': str(self.hitAccuracy),
                    'ruemmler': str(self.ruemmler_stability_criterion),
                    'nyquist_criteria': str(self.open_loop_stability)
                },
                'General settings': {
                    # 'recalcCharEq': str(self.recalcCharEq),
                    'Frequency range': str(self.freq_range),
                    'Frequency step': str(self.df),
                    'Time lag range oxidizer': str(self.tau_Ox),
                    'Time lag step oxidizer': str(self.dtau_Ox),
                    'Time lag range fuel': str(self.tau_Fu),
                    'Time lag step fuel': str(self.dtau_Fu),
                    'Interaction index': str(self.N),
                    'Interaction index step': str(self.dN)
                },
                'Chamber geometry': {
                    'Chamber diameter': str(self.d_cc),
                    'Chamber length': str(self.l_cc),
                    'Contraction ratio': str(self.epsilon_cc)
                },
                'Injector geometry': {
                    'Config File': str(self.injector_configuration)
                },
                'Load point definition': {
                    'Fuel': str(self.fuel),
                    'Chamber pressure': str(self.p_cc),
                    'ROF primary injection': str(self.rof_cc),
                    'LOX mass flow': str(self.mDot_O2),
                    'Fuel mass flow': str(self.mDot_F),
                    'Secondary fuel mass flow': str(self.mDot_SF),
                    'Window cooling fluid': str(self.fluid_wc),
                    'Window cooling mass flow': str(self.mDot_WC),
                    'LOX dome pressure': str(self.p_io),
                    'Fuel dome pressure': str(self.p_if),
                    'Secondary fuel dome pressure': str(self.p_sf),
                    'Window cooling dome pressure': str(self.p_wc),
                    'LOX injection temperature': str(self.T_io),
                    'Fuel injection temperature': str(self.T_if),
                    'Window cooling injection temperature': str(self.T_wc),
                    'Combustion efficiency': str(self.eta_cc)
                },
                'Frequencies in tests': {
                    'Main frequency': str(self.mainFreq),
                    'Frequencies in tests': str(self.freq_tests)
                }
            },
                default_flow_style=False, indent=4, sort_keys=False)

            # Write to yaml_file while doing some pretty formatting by removing all ' characters
            # (PyYaml creates a lot of them)
            yaml_file.write(yaml_content.replace("'", "").replace('nan', '-').replace('[]', '-'))

        with open(os.path.join(self.result_path, INPUT_CONFIGURATIONS_DIRECTORY, CONFIG_FILE_NAME + '.csv'), 'w',
                  newline='') as csv_file:
            header = ['# Identifier', self.identifier]
            writer = csv.DictWriter(csv_file, fieldnames=header)
            writer.writeheader()
            for attr, value in self.__dict__.items():
                if not any(k == attr for k in ['identifier', 'result_path', 'injector', 'combustion_chamber']):
                    writer.writerow({'# Identifier': attr, self.identifier: value})

    def reload_injector(self):
        return injector_configuration(self, injector_filename=self.injector_configuration, fluid_fuel=self.fuel,
                                      fluid_window_cooling=self.fluid_wc,
                                      T_io=self.T_io, T_if=self.T_if, T_sf=self.T_sf, T_wc=self.T_wc,
                                      mDot_O2=self.mDot_O2, mDot_F=self.mDot_F, mDot_SF=self.mDot_SF,
                                      mDot_WC=self.mDot_WC)


def handle_file_deleting_error(func, path, exc_info):
    """ Windows is bad.

        Sometimes, Windows thinks it's a good idea to create files with read-only permissions. In this case, python
        will fail to delete them (they're read-only). However, shutil.rmtree can try to handle this error by changing
        the permission manually, and then trying to delete the files again.

    """
    if not os.access(path, os.W_OK):
        # Change the permission to write
        # print("Hey")
        os.chmod(path, stat.S_IWUSR)
        # Call the function again to delete file finally
        func(path)
