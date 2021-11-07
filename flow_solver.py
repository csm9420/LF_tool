"""
Created on Thu Jul 14 08:49:45 2021

@author: spah_mr
Based on the bachelor thesis of Malte Gradetzke. View flow_solver.m.
The initial assumed values for the iteration are changed.

# Defining the symbols used here
    #   p           local static pressure [Pa]
    #   p_total     total pressure [Pa]
    #   dp          local static pressure loss [Pa]
    #   dp_total    local total pressure loss [Pa]
    #   T           local static temperature [K]
    #   T_total     total temperature [K]
    #   rho         local density [kg/m^3]
    #   rho_total   total density [kg/m^3]
    #   c           local speed of sound [m/s]
    #   h           local enthalpy [J/kg]
    #   h_total     total enthalpy [J/kg]
    #   u           internal energy [J/kg]
    #   mdot        mass flow [kg/s]
    #   zeta        loss coefficient for the total pressure loss over an element [-]
    #   kappa       adiabatic index [-]
    #   D           diameter [m]
    #   Da          outer diameter [m]
    #   Di          inner diameter [m]
    #   Dh          hydraulic diameter [m]
    #   eta         dynamic viscosity [Pa s]
    #   A           Cross section area [m^2]
    #   V           Volume [m^3]

"""

from Refprop import getRefPropValue
from math import pow, sqrt, log10, pi
import pandas as pd
from Injector import new_manifold, new_orifice, new_line
import os.path
from sys import exit
from tqdm import tqdm
import matplotlib.pyplot as plt
from numpy import isnan

# Gas constant, see https://en.wikipedia.org/wiki/Gas_constant
GAS_CONSTANT = 8.31446261815324  # [J/(K mol)]
# Molar mass oxygen, see https://pubchem.ncbi.nlm.nih.gov/compound/oxygen
MOLAR_MASS_OXYGEN = 0.031999  # [kg/mol]
# Molar mass hydrogen, see https://pubchem.ncbi.nlm.nih.gov/compound/783
MOLAR_MASS_HYDROGEN = 0.00216  # [kg/mol]

# Convergence limit for the total pressure drop iteration
DP_TOTAL_ITERATION_ACCURACY = 1e-10
# Convergence limit for the local static density iteration
RHO_ITERATION_ACCURACY = 1e-10
# RHO_ITERATION_ACCURACY = 1e-1
# Maximum allowed iterations
ITERATION_MAX = 200
# Number of nodes for the total pressure loss calculation within the tube of each element per mm length
AMOUNT_OF_NODES_PER_MM = 4


# Defining the ranges (= the rows) in the excel injector configuration file
# Do not change this, unless the excel file is changed
RANGE_LOX_INJECTOR = range(14, 24)
RANGE_FUEL_INJECTOR = range(38, 48)
RANGE_LOX_FEED = range(2, 12)
RANGE_FUEL_FEED = range(26, 36)
RANGE_SECONDARY_FUEL_INJECTOR = range(50, 60)
RANGE_WINDOW_COOLING_INJECTOR = range(62, 72)


class FlowSolverInjector:

    def __init__(self, path, settings_file, sheet, mDot_LOX, mDot_fuel, p_io, p_if, T_io, T_if, p_feed_o1, p_feed_f1,
                 p_feed_o2, p_feed_f2, T_feed_f1, T_feed_o1, T_feed_f2, T_feed_o2, p_cc, fuel_type, p_sf, T_sf, p_wc, T_wc,
                 mDot_sf, mDot_wc, fluid_wc, use_feed_system):
        """ Flow solver

        :param path: path to the result file
        :param settings_file: Excel settings file
        :param sheet: sheet within settings_file of the injector
        :param mDot_LOX: oxygen mass flow [kg/s]
        :param mDot_fuel: fuel mass flow [kg/s]
        :param p_io: pressure oxygen injector manifold [Pa]
        :param p_if: pressure fuel injector manifold [Pa]
        :param T_io: temperature oxygen injector manifold [K]
        :param T_if: temperature fuel injector manifold [K]
        :param p_feed_o1: Feed system oxygen manifold pressure [Pa]
        :param p_feed_f1: Feed system fuel manifold pressure [Pa]
        :param T_feed_o1: Feed system oxygen manifold pressure [K]
        :param T_feed_f1: Feed system fuel manifold pressure [K]
        :param p_feed_o2: Feed system oxygen manifold pressure [Pa]
        :param p_feed_f2: Feed system fuel manifold pressure [Pa]
        :param T_feed_o2: Feed system oxygen manifold pressure [K]
        :param T_feed_f2: Feed system fuel manifold pressure [K]
        :param p_cc: pressure combustion chamber [Pa]
        :param fuel_type: type of fuel
        :param p_sf: Pressure of the secondary fuel manifold [Pa]
        :param T_sf: Temperature of the secondary fuel manifold [K]
        :param p_wc: Pressure of the window cooling manifold [Pa]
        :param T_wc: Temperature of the window cooling manifold [K]
        :param mDot_sf: Mass flow secondary fuel [kg/s]
        :param mDot_wc: Mass flow window cooling [kg/s]
        :param fluid_wc: Defines the fluid used for the window cooling
        :param use_feed_system: Determines, if the feed system should be used. If False, use only the injector part
        :return: tuple containing lists of the injector elements for the LF tool

        :type path: os.path
        :type settings_file: str
        :type sheet: str
        :type mDot_LOX: float
        :type mDot_fuel: float
        :type p_io: float
        :type p_if: float
        :type T_io: float
        :type T_if: float
        :type p_cc: float
        :type p_feed_o1: float
        :type T_feed_o1: float
        :type T_feed_f1: float
        :type p_feed_f1: float
        :type p_feed_o2: float
        :type T_feed_o2: float
        :type T_feed_f2: float
        :type p_feed_f2: float
        :type fuel_type: str
        :type p_sf: float
        :type T_sf: float
        :type p_wc: float
        :type T_wc: float
        :type mDot_sf: float
        :type mDot_wc: float
        :type fluid_wc: str
        :type use_feed_system: bool
        """

        # Define the initial attributes
        self.LOX = None
        self.Fuel = None

        # Set the fuel type correctly (for Coolprop)
        # TODO: More fuel types needed?
        if fuel_type == 'H2':
            fuel_type = 'Hydrogen'

        # Set the window cooling type correctly (for Coolprop)
        if fluid_wc == 'H2':
            fluid_wc = 'Hydrogen'

        # Define the settings file
        file = os.path.join(os.getcwd(), settings_file)

        # Create the FlowElements based on the settings file
        LOX_injector = create_flow_elements(file, sheet, path, 'LOX_injector', RANGE_LOX_INJECTOR, mDot_LOX)
        Fuel_injector = create_flow_elements(file, sheet, path, 'Fuel_injector', RANGE_FUEL_INJECTOR, mDot_fuel)
        Secondary_fuel_injector = create_flow_elements(file, sheet, path, 'Secondary_Fuel_Injector', RANGE_SECONDARY_FUEL_INJECTOR, mDot_sf)
        Window_cooling_injector = create_flow_elements(file, sheet, path, 'Window_Cooling_Injector', RANGE_WINDOW_COOLING_INJECTOR, mDot_wc)

        # Determines if the feed system should be used, independent if it's specified or not.
        if use_feed_system:
            LOX_feed = create_flow_elements(file, sheet, path, 'LOX_feed', RANGE_LOX_FEED, mDot_LOX)
            Fuel_feed = create_flow_elements(file, sheet, path, 'Fuel_feed', RANGE_FUEL_FEED, mDot_fuel)
        else:
            LOX_feed = []
            Fuel_feed = []

        # Create the injector elements by calculating the FlowElements
        LOX_injector = calculate_elements(path, 'LOX_injector', LOX_injector, p_io, T_io, p_cc, 'Oxygen')
        Fuel_injector = calculate_elements(path, 'Fuel_injector', Fuel_injector, p_if, T_if, p_cc, fuel_type)

        # Create the feed system (injector) elements by calculating the FlowElements
        # This calculation is only done if there are any elements defined in the settings file. If not, this will stay
        # an empty list
        # TODO: Find a source for this assumption.
        #   This uses p_io and p_if as p_end for the feed system parts. There's no source for that found.
        #   However, if the injector part (the actual injector) uses p_cc as p_end (-> the pressure of the manifold in
        #   which the system is injecting), then it should be p_io and p_if for the feed system parts.
        if LOX_feed:
            LOX_feed = calculate_elements_feed(path, 'LOX_feed', LOX_feed, p_feed_o1, p_feed_o2, T_feed_o1, T_feed_o2, p_io, p_cc, 'Oxygen')
        if Fuel_feed:
            Fuel_feed = calculate_elements_feed(path, 'Fuel_feed', Fuel_feed, p_feed_f1, p_feed_f2, T_feed_f1, T_feed_f2, p_if, p_cc, fuel_type)
        if Secondary_fuel_injector:
            Secondary_fuel_injector = calculate_elements(path, 'Secondary_fuel_injector', Secondary_fuel_injector,
                                                         p_sf, T_sf, p_cc, fuel_type)
        if Window_cooling_injector:
            Window_cooling_injector = calculate_elements(path, 'Window_cooling_injector', Window_cooling_injector,
                                                         p_wc, T_wc, p_cc, fluid_wc)

        # Connect feed system and injector
        LOX_system = LOX_feed + LOX_injector
        Fuel_system = Fuel_feed + Fuel_injector

        self.LOX = LOX(LOX_system)
        self.Fuel = Fuel(Fuel_system, fuel_type)
        if Secondary_fuel_injector:
            self.SecFuel = SecFuel(Secondary_fuel_injector, fuel_type)
        if Window_cooling_injector:
            self.WinCool = WinCool(Window_cooling_injector, fluid_wc)

    # Define Function as needed to create the transfer function from the original injector definition
    def getLengthCas(self, fluid):
        accumLength = 0
        if fluid == "Oxygen":
            for i in self.LOX.elements:
                if i.type == "Line":
                    accumLength += i.length
        elif fluid == "Fuel":
            for i in self.Fuel.elements:
                if i.type == "Line":
                    accumLength += i.length
        return accumLength


class LOX:

    def __init__(self, elements):
        self.num_elements = len(elements)
        self.elements = elements
        self.type = 'LOX'


class Fuel:

    def __init__(self, elements, fluid):
        self.num_elements = len(elements)
        self.elements = elements
        self.type = fluid


class SecFuel:

     def __init__(self, elements, fluid):
         self.num_elements = len(elements)
         self.elements = elements
         self.type = fluid


class WinCool:

    def __init__(self, elements, fluid):
        self.num_elements = len(elements)
        self.elements = elements
        self.type = fluid


class FlowElement:

    def __init__(self, data, mDot):
        """ Initialize a single FlowElement

        :param data: a pandas Series containing a the information about the single FlowElement
        :type data: ps.Series
        :rtype: FlowElement

        """
        self.Da = data.loc['outer diameter'] * 1e-3
        self.Di = data.loc['inner diameter'] * 1e-3
        self.Dh = data.loc['hydraulic diameter'] * 1e-3
        # self.Dh = self.Da - self.Di
        self.A = data.loc['cross section area'] * 1e-6
        # self.A = (pow(self.Da, 2) - pow(self.Di, 2)) * pi / 4
        self.L = data.loc['length'] * 1e-3
        self.V = data.loc['volume'] * 1e-9
        self.k = data.loc['surface roughness'] * 1e-3
        self.zeta = data.loc['inlet loss coefficient']
        self.node_length = None
        self.nodes = None
        self.number_of_elements = data.loc['number of elements']
        self.mDot = mDot / self.number_of_elements
        self.type = data.loc['type']
        # self.h_total = None
        # self.T_total = None
        # self.p_total = None
        # self.rho_total = None
        # self.h = None
        # self.T = None
        # self.p = None
        # self.rho = None
        # self.dp_total = None
        # self.dp = None
        # self.v = None
        # self.c = None
        # self.Ma = None
        # self.Re = None
        # self.eta = None
        # self.Lamda = None
        self.flow = None

    def create_nodes(self, offset):
        """ Create all nodes for the current FlowElement, where the value of a node refers to a absolute x position.
            Per definition, the end of the injector is defined as x = 0 [m].

                      element_{i-1}                |            element_{i}
            =======================================|=======================================
               0    1    2   3   4   5   6   7    8|   0    1    2   3   4   5   6   7    8
            =======================================|=======================================
                                                       |                                     x=0


        :param offset: offset of the end of the FlowElement from x = 0 [m]
        :type offset: float
        """

        length = self.L * 1e3  # Length of the current element [mm]
        number_of_nodes = round(length * AMOUNT_OF_NODES_PER_MM)  # number of nodes for the current element
        self.node_length = self.L / number_of_nodes  # Distance between two nodes [m]

        nodes = []
        for i in reversed(range(number_of_nodes)):
            nodes.append(0 - offset - i * self.node_length)

        # Nodes for the current element, where a node refers to a x position [m]
        self.nodes = nodes

def create_flow_elements(file, sheet, path, name, location, mDot):
    """ Creates a list of FlowElements, ready for flow parameter calculation.
        The list is required to start with a manifold, followed by pipe elements.

        :param file: path to the excel file where the information is stored
        :param location: Location/Range of cells in the excel file
        :param path: path to store some files
        :param name: Name of the current element system
        :param mDot: Total mass flow through the fluid system [kg/s]
        :return: A list of FlowElements containing all input parameters
        :type file: os.path
        :rtype: list[FlowElement]
    """

    # Read the range from the excel file, rename the variables to the identifiers, drop them, drop all empty columns
    # data is then a pd.DataFrame containing information about all flow elements
    data = pd.read_excel(file, sheet_name=sheet, engine='openpyxl')
    data = data.iloc[location].rename(index=lambda k: data.loc[k]['# Identifier']).drop('# Identifier', axis='columns') \
        .drop('# Unit', axis='columns').drop('# Symbol', axis='columns').dropna(axis='columns', how='all')

    elements = [FlowElement(data[j], mDot) for j in data]
    # Assign the nodes with the respective offset for each element, ignoring the first element (manifold)
    offset = 0
    if (len(elements) > 0 and elements[0].type != "Manifold"):
        for i in reversed(range(len(elements))):
            elements[i].create_nodes(offset)
            offset += elements[i].L
    else:
        for i in reversed(range(len(elements) - 1)):
            elements[i + 1].create_nodes(offset)
            offset += elements[i + 1].L

    # Calculate zeta values, if needed
    for i, element in enumerate(elements):
        if element.zeta == 'auto':
            element.zeta = predict_zeta(elements[i - 1].Dh, elements[i].Dh)

    elements_to_file(path, name, elements)

    return elements


def predict_zeta(Dh_prev, Dh_current):
    """ Predicts zeta values, if not provided.

    :param Dh_prev: Hydraulic diameter of the previous element
    :param Dh_current: Hydraulic diameter of the current element
    :return: Zeta value for the inlet of the current element
    :type Dh_current: float
    :type Dh_prev: float
    :rtype: float
    """

    if isnan(Dh_prev):
        # The previous element is a manifold which has no hydraulic diameter
        # Not able to predict a value
        # zeta = 0.5
        print('Not able to predict a zeta value if the previous element is manifold because a hydraulic diameter '
              'is needed.\nProvide a value like zeta = 0.5 manually. EXIT.')
        exit()
    elif Dh_prev < Dh_current:
        zeta = (1 - (Dh_prev / Dh_current) ** 2) ** 2
    else:
        zeta = 0.5 * (1 - (Dh_current / Dh_prev) ** 2)

    return zeta

def calculate_elements_feed(path, name, elements, p1, p2, T1, T2, p_end, p_cc, fluid):
    """ Calculates flow mechanics parameters for feeding system utilizing feed system information
    acquired from W. Armbruster - S. Choi
    :param path: Path to a folder to store the results
    :param elements: List of FlowElements, starting with a a manifold
    :param p1, p2 : pressure at the manifold and line connection
    :param T1, T2 : temperature at the manifold and line connection
    :param p_end: Static pressure [Pa] after the last element
    """

    # Initialization of injector elements (Modified elements matrix including orifices between line elements)
    injector_elements = []
    # Create an object to monitor the progress; this is just for display and monitoring purposes
    monitor = tqdm(total=len(elements),
                   desc='Calculate flow for {}'.format(name), leave=True,
                   bar_format='{desc:<65}| {bar:25}| {n_fmt}/{total_fmt}')

    elements[0].flow = pd.DataFrame(data=
                                    {'x_nodes': 'manifold',
                                     'h_total': getRefPropValue('H', 'T', T1, 'P', p1, fluid),
                                     'u': getRefPropValue('U', 'T', T1, 'P', p1, fluid),
                                     'T_total': T1,
                                     'T': T1,
                                     'p_total': p1,
                                     'p': p1,
                                     'dp_total': 0,
                                     'rho_total': getRefPropValue('D', 'T', T1, 'P', p1, fluid),
                                     'rho': getRefPropValue('D', 'T', T1, 'P', p1, fluid),
                                     'c': getRefPropValue('A', 'T', T1, 'P', p1, fluid),
                                     'eta': getRefPropValue('V', 'T', T1, 'P', p1, fluid),
                                     'Lambda': None,
                                     'v': 0,
                                     'Re': 0,
                                     'Ma': 0}, index=[0])
    # Setting basic variables and acoustic parameters
    injector_elements.append(new_manifold(elements[0].number_of_elements, elements[0].mDot, elements[0].V,
                                          elements[0].flow['c'].iloc[-1],
                                          elements[0].flow['rho'].iloc[-1], p_cc))
    injector_elements[0].p = p1

    # Update monitoring, just for display purposes
    monitor.update()

    # Iterating over all elements but skipping the first one because the parameters are already determined
    for i in range(1, len(elements)):
        if (i == 1):
            # Pressure drop at the orifice to calculate pressure at the start of the 1st line
            dp_total_inlet = total_pressure_loss_inlet(path=path, mdot_1=elements[i].mDot,
                                                       p_total_0=elements[i - 1].flow['p_total'].iloc[-1],
                                                       rho_total_0=elements[i - 1].flow['rho_total'].iloc[-1],
                                                       rho_0=elements[i - 1].flow['rho'].iloc[-1],
                                                       T_total_0=elements[i - 1].flow['T_total'].iloc[-1],
                                                       h_total_0=elements[i - 1].flow['h_total'].iloc[-1],
                                                       zeta_1=elements[i].zeta,
                                                       A_1=elements[i].A,
                                                       fluid=fluid)
            # Calculate the total pressure of the first tube element / at the start of the tube
            p_total_tube_init = elements[i - 1].flow['p_total'].iloc[-1] - dp_total_inlet
            # Element tube
            elements[i].flow = tube_conditions(path=path,
                                   mdot=elements[i].mDot,
                                   pt_1=p_total_tube_init,
                                   T_total_0=elements[i - 1].flow['T_total'].iloc[-1],
                                   h_total_0=elements[i - 1].flow['h_total'].iloc[-1],
                                   A=elements[i].A,
                                   Di=elements[i].Di,
                                   Dh=elements[i].Dh,
                                   node_length=elements[i].node_length,
                                   nodes=elements[i].nodes,
                                   k=elements[i].k,
                                   fluid=fluid)

        else:
            # Pressure and temperature at the start of the 2nd line is measured as P2 and T2
            elements[i].flow = tube_conditions(path=path,
                                   mdot=elements[i].mDot,
                                   pt_1=p2,
                                   T_total_0=T2,
                                   h_total_0=elements[i - 1].flow['h_total'].iloc[-1],
                                   A=elements[i].A,
                                   Di=elements[i].Di,
                                   Dh=elements[i].Dh,
                                   node_length=elements[i].node_length,
                                   nodes=elements[i].nodes,
                                   k=elements[i].k,
                                   fluid=fluid)

        # Add the orifice element to injector_elements
        if (i == 1):
            dp_orifice = p1 - p_total_tube_init
        else:
            dp_orifice = p_total_tube_init - p2
        # Orifice element addition to return value - injector_elements
        injector_elements.append(new_orifice(elements[i].number_of_elements, elements[i].mDot, dp_orifice, p_cc))
        # Pressure value at the orifice element
        injector_elements[i*2-1].p = injector_elements[i*2-1-1].p - injector_elements[i*2-1].dp
        # After the orifice, add a line element
        # Because v, c are not perfectly constant along the line, use mean values
        tube_mean_v = elements[i].flow['v'].mean()
        tube_mean_c = elements[i].flow['c'].mean()
        # Line element addition behind the orifice (next line)
        injector_elements.append(new_line(elements[i].number_of_elements, elements[i].mDot, elements[i].A, tube_mean_v,
                                          tube_mean_c, elements[i].L, p_cc))
        # Pressure value at the line element
        injector_elements[i*2].p = injector_elements[i*2-1].p
        # Update monitoring, just for display purposes
        monitor.update()

    # Collect all results in a Dataframe and save in a .csv file
    flow = pd.DataFrame(columns=['x_nodes', 'h_total', 'u', 'T_total', 'T', 'p_total', 'p', 'dp_total',
                                 'rho_total', 'rho', 'c', 'eta', 'Lambda', 'v', 'Re', 'Ma'])
    for element in elements:
        flow = flow.append(element.flow, ignore_index=True)
    flow.to_csv(os.path.join(path, '{}_flow_solver.csv'.format(name)), index=False)

    # Plot
    plot_pressure(os.path.join(path, '{}_pressure.png'.format(name)), name, flow['x_nodes'],
                  flow['p'], flow['p_total'], p1, p_end)
    plot_velocity(os.path.join(path, '{}_velocity.png'.format(name)), name, flow['x_nodes'], flow['v'])
    plot_speed_of_sound(os.path.join(path, '{}_speed_of_sound.png'.format(name)), name, flow['x_nodes'], flow['c'])
    plot_temperature(os.path.join(path, '{}_temperature.png'.format(name)), name, flow['x_nodes'], flow['T'],
                     flow['T_total'], T1)
    plot_density(os.path.join(path, '{}_density.png'.format(name)), name, flow['x_nodes'], flow['rho'],
                 flow['rho_total'])
    # plot_geometry(os.path.join(path, '{}_geometry.png'.format(name)), name, flow['x_nodes'], flow['Da'], flow['Di'])

    # Close monitoring, just for display purposes
    monitor.close()

    # Add a last element to end with an orifice, whereby the pressure drop joins the calculated pressure and the
    # end pressure of the system, here feed system
    dp_last_orifice = p2 - p_end
    # Velocity at the orifice
    v_last_orifice = elements[-1].flow['v'].iloc[-1]
    injector_elements.append(new_orifice(number_of_elements=elements[-1].number_of_elements, mDot=elements[-1].mDot,
                                         dp=dp_last_orifice, p_end=p_cc, v=v_last_orifice))
    # Pressure at the last orifice between the 2nd feed line and the injection manifold
    injector_elements[-1].p = p_end
    # Return the injector elements for the LF Tool
    return injector_elements

def calculate_elements(path, name, elements, p_total_0, T_total_0, p_end, fluid):
    """ Calculates all flow mechanics parameters for a list of FlowElements, whereby the first element has to be a
    manifold.

    :param path: Path to a folder to store the results
    :param elements: List of FlowElements, starting with a a manifold
    :param p_total_0: Total pressure [Pa] of the first element
    :param T_total_0: Total temperature [K] of the first element
    :param p_end: Static pressure [Pa] after the last element
    :param fluid: Type of the current fluid
    :type elements: list[FlowElement]
    :type p_total_0: float
    :type T_total_0: float
    :type p_end: float
    :type fluid: str
    :rtype: list[Injector.Element]
    """
    #
    # Assumptions
    # - The first element is a manifold.
    # - The feed system is adiabatic.
    # - The system is in a steady state.
    #
    #
    #
    # A single element is separated in the element inlet and a following tube. There's a total pressure loss over the
    # inlet, and an additional total pressure loss over the tube of the element.

    # Creating a list of injector elements which will be returned
    injector_elements = []

    # Create an object to monitor the progress; this is just for display and monitoring purposes
    monitor = tqdm(total=len(elements),
                   desc='Calculate flow for {}'.format(name), leave=True,
                   bar_format='{desc:<65}| {bar:25}| {n_fmt}/{total_fmt}')
    if(elements[0].type == "Manifold"):
        elements[0].flow = pd.DataFrame(data=
                                        {'x_nodes': 'manifold',
                                         'h_total': getRefPropValue('H', 'T', T_total_0, 'P', p_total_0, fluid),
                                         'u': getRefPropValue('U', 'T', T_total_0, 'P', p_total_0, fluid),
                                         'T_total': T_total_0,
                                         'T': T_total_0,
                                         'p_total': p_total_0,
                                         'p': p_total_0,
                                         'dp_total': 0,
                                         'rho_total': getRefPropValue('D', 'T', T_total_0, 'P', p_total_0, fluid),
                                         'rho': getRefPropValue('D', 'T', T_total_0, 'P', p_total_0, fluid),
                                         'c': getRefPropValue('A', 'T', T_total_0, 'P', p_total_0, fluid),
                                         'eta': getRefPropValue('V', 'T', T_total_0, 'P', p_total_0, fluid),
                                         'Lambda': None,
                                         'v': 0,
                                         'Re': 0,
                                         'Ma': 0}, index=[0])
        injector_elements.append(new_manifold(elements[0].number_of_elements, elements[0].mDot, elements[0].V,
                                              elements[0].flow['c'].iloc[-1],
                                              elements[0].flow['rho'].iloc[-1], p_end))
        injector_elements[0].p = p_total_0
    else:
        elements[0].flow = tube_conditions(path=path,
                               mdot=elements[0].mDot,
                               pt_1=p_total_0,
                               T_total_0=T_total_0,
                               h_total_0=getRefPropValue('H', 'T', T_total_0, 'P', p_total_0, fluid),
                               A=elements[0].A,
                               Di=elements[0].Di,
                               Dh=elements[0].Dh,
                               node_length=elements[0].node_length,
                               nodes=elements[0].nodes,
                               k=elements[0].k,
                               fluid=fluid)
        tube_mean_v = elements[0].flow['v'].mean()
        tube_mean_c = elements[0].flow['c'].mean()
        injector_elements.append(new_line(elements[0].number_of_elements, elements[0].mDot, elements[0].A, tube_mean_v,
                                          tube_mean_c, elements[0].L, p_end))
        injector_elements[0].p = p_total_0

    # Update monitoring, just for display purposes
    monitor.update()

    # Iterating over all elements but skipping the first one because the parameters are already determined
    for i in range(1, len(elements)):
        # Element Inlet
        dp_total_inlet = total_pressure_loss_inlet(path=path, mdot_1=elements[i].mDot,
                                                   p_total_0=elements[i - 1].flow['p_total'].iloc[-1],
                                                   rho_total_0=elements[i - 1].flow['rho_total'].iloc[-1],
                                                   rho_0=elements[i - 1].flow['rho'].iloc[-1],
                                                   T_total_0=elements[i - 1].flow['T_total'].iloc[-1],
                                                   h_total_0=elements[i - 1].flow['h_total'].iloc[-1],
                                                   zeta_1=elements[i].zeta,
                                                   A_1=elements[i].A,
                                                   fluid=fluid)

        # Element tube
        # Calculate the total pressure of the first tube element / at the start of the tube
        p_total_tube_init = elements[i - 1].flow['p_total'].iloc[-1] - dp_total_inlet

        # Calculate tube conditions, here, tube is a data frame
        tube = tube_conditions(path=path,
                               mdot=elements[i].mDot,
                               pt_1=p_total_tube_init,
                               T_total_0=elements[i - 1].flow['T_total'].iloc[-1],
                               h_total_0=elements[i - 1].flow['h_total'].iloc[-1],
                               A=elements[i].A,
                               Di=elements[i].Di,
                               Dh=elements[i].Dh,
                               node_length=elements[i].node_length,
                               nodes=elements[i].nodes,
                               k=elements[i].k,
                               fluid=fluid)

        elements[i].flow = tube

        # Add the orifice element to injector_elements
        # Definition: If the local static pressure behind the orifice is lower than before, than dp_orifice > 0
        # TODO: It is not totally clear whether the static pressure or the total pressure should be used here
        # dp_orifice = elements[i - 1].flow['p'].iloc[-1] - elements[i].flow['p'].iloc[0]
        dp_orifice = elements[i - 1].flow['p_total'].iloc[-1] - elements[i].flow['p_total'].iloc[0]
        injector_elements.append(new_orifice(elements[i].number_of_elements, elements[i].mDot, dp_orifice, p_end))
        injector_elements[i*2-1].p = injector_elements[i*2-1-1].p - injector_elements[i*2-1].dp
        # After the orifice, add a line element
        # Because v, c are not perfectly constant along the line, use mean values
        tube_mean_v = tube['v'].mean()
        tube_mean_c = tube['c'].mean()
        injector_elements.append(new_line(elements[i].number_of_elements, elements[i].mDot, elements[i].A, tube_mean_v,
                                          tube_mean_c, elements[i].L, p_end))
        injector_elements[i * 2].p = injector_elements[i * 2 - 1].p

        # Update monitoring, just for display purposes
        monitor.update()

    # Collect all results in a Dataframe and save in a .csv file
    flow = pd.DataFrame(columns=['x_nodes', 'h_total', 'u', 'T_total', 'T', 'p_total', 'p', 'dp_total',
                                 'rho_total', 'rho', 'c', 'eta', 'Lambda', 'v', 'Re', 'Ma'])
    for element in elements:
        flow = flow.append(element.flow, ignore_index=True)
    flow.to_csv(os.path.join(path, '{}_flow_solver.csv'.format(name)), index=False)

    # Plot
    plot_pressure(os.path.join(path, '{}_pressure.png'.format(name)), name, flow['x_nodes'],
                  flow['p'], flow['p_total'], p_total_0, p_end)
    plot_velocity(os.path.join(path, '{}_velocity.png'.format(name)), name, flow['x_nodes'], flow['v'])
    plot_speed_of_sound(os.path.join(path, '{}_speed_of_sound.png'.format(name)), name, flow['x_nodes'], flow['c'])
    plot_temperature(os.path.join(path, '{}_temperature.png'.format(name)), name, flow['x_nodes'], flow['T'],
                     flow['T_total'], T_total_0)
    plot_density(os.path.join(path, '{}_density.png'.format(name)), name, flow['x_nodes'], flow['rho'],
                 flow['rho_total'])
    # plot_geometry(os.path.join(path, '{}_geometry.png'.format(name)), name, flow['x_nodes'], flow['Da'], flow['Di'])

    # Close monitoring, just for display purposes
    monitor.close()

    # Add a last element to end with an orifice, whereby the pressure drop joins the calculated pressure and the
    # end pressure of the system
    # Pressure drop to the chamber
    # TODO: It is not totally clear whether the static pressure or the total pressure should be used here
    #   If the static pressure is used, p_end can be used because the measured pressure is the static pressure
    #   However, if the total pressure is used, the total pressure of the chamber has to be used
    #   (which is p_end_total = p_end + p_end_dynamic)
    # dp_last_orifice = elements[-1].flow['p'].iloc[-1] - p_end
    dp_last_orifice = elements[-1].flow['p_total'].iloc[-1] - p_end
    # Velocity at the orifice
    v_last_orifice = elements[-1].flow['v'].iloc[-1]
    injector_elements.append(new_orifice(number_of_elements=elements[-1].number_of_elements, mDot=elements[-1].mDot,
                                         dp=dp_last_orifice, p_end=p_end, v=v_last_orifice))
    injector_elements[-1].p = p_end

    # Return the injector elements for the LF Tool
    return injector_elements


def total_pressure_loss_inlet(path, mdot_1, p_total_0, rho_total_0, rho_0, T_total_0, h_total_0, zeta_1, A_1, fluid):
    """ Calculates the total pressure loss over the current element (1).

    The inlet of the current element is 1, the element before is 0, the end of the element is 2.

    ==========
              |==================
      0       1       2
              |==================
    ==========

    :param path: path to store fome files
    :param mdot_1: mass flow [kg/s]
    :param p_total_0: total pressure before the inlet [Pa]
    :param rho_total_0: total density before the inlet [kg/m^3]
    :param rho_0: static density before the inlet [kg/m^3]
    :param T_total_0: total temperature before the inlet [K]
    :param h_total_0: total enthalpy before the inlet [J/kg]
    :param zeta_1: loss coefficient of the inlet [-]
    :param A_1: cross section area of the inlet [m^2]
    :param fluid: specifies the current fluid
    :return: total pressure drop over the inlet [Pa]
    :type path: os.path
    :type mdot_1: float
    :type p_total_0: float
    :type rho_total_0: float
    :type rho_0: float
    :type T_total_0: float
    :type h_total_0: float
    :type zeta_1: float
    :type A_1: float
    :type fluid: str
    :rtype: float
    """

    # The local static density is needed, but not know. Assume a start value.
    # TODO: Use a better estimation (how?)
    # rho = rho_total_0  # Density [kg/m^3]
    rho = rho_0  # Density [kg/m^3]

    # Assuming a steady state, the velocity can be determined
    v = mdot_1 / (rho * A_1)
    # Using the values based on the guessed density, the initial guess for the total pressure drop between element 0
    # and element 1 is calculated.
    dp_total_iter_init = zeta_1 * (rho / 2) * pow(v, 2)

    # Set the total pressure loss delta to an initial value. This should converge during the iteration.
    delta = DP_TOTAL_ITERATION_ACCURACY * 1.1
    # Set the iteration counter to its initial values. This will prevent an endless loop.
    iter_count = 0

    while abs(delta) > DP_TOTAL_ITERATION_ACCURACY and iter_count < ITERATION_MAX:
        # Calculate the total pressure behind the current element
        p_total_2 = p_total_0 - dp_total_iter_init

        p_init_adibatic_tube = p_total_2 - rho * pow(v, 2) / 2

        # Calculate the density behind the element assuming an adiabatic tube element
        _, _, rho, v, _, _ = adiabatic_tube(path, mdot_1, p_total_2, p_init_adibatic_tube, rho_total_0, rho, T_total_0,
                                            A_1, fluid)

        # Using the newly acquired values for the density and the velocity, re-calculate the total pressure loss
        dp_total_iter_new = zeta_1 * (rho / 2) * pow(v, 2)

        # Determine the convergence, update the initial density value for the next iteration step, updating the
        # iteration counter
        delta = dp_total_iter_new - dp_total_iter_init
        dp_total_iter_init = dp_total_iter_new
        iter_count += 1

    # Print an iteration waring, if necessary
    iteration_max_warning(iter_count, 'total_pressure_loss_inlet')

    # Use the final parameters after the iteration
    dp_total = dp_total_iter_new

    return dp_total

def tube_conditions(path, mdot, pt_1, T_total_0, h_total_0, A, Di, Dh, k, nodes, node_length, fluid):
    """

    :param path: path to store some files
    :param mdot: mass flow [kg/s]
    :param pt_1: total pressure after the total inlet pressure loss [Pa]
    :param T_total_0: total temperature before the inlet [K]
    :param h_total_0: total enthalpy before the inlet [J/kg]
    :param A: cross section area of the tube/inlet [m^2]
    :param Di: inner diameter [m]
    :param Dh: hydraulic diameter (dh = da - di) [m^2]
    :param k: surface roughness [m]
    :param nodes: list of all x-positions of nodes
    :param node_length: length of a node [m]
    :param fluid: specifies the current fluid
    :return: DataFrame containing all flow parameters

    :type path: os.path
    :type mdot: float
    :type pt_1: float
    :type T_total_0: float
    :type h_total_0: float
    :type A: float
    :type Di: float
    :type Dh: float
    :type k: float
    :type nodes: list[float]
    :type node_length: float
    :type fluid: str
    :rtype: pd.DataFrame
    """

    # Create lists to collect all values
    h_total, u, T_total, p_total, T, p_total, p, dp_total, rho_total, rho, c, eta, Lambda, v, Re, Ma = \
        [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

    # Set the initial values for the total parameters
    h_total.append(h_total_0)
    T_total.append(T_total_0)
    p_total.append(pt_1)

    # Iterating over all nodes. Here, x_node is the x position of the node, and nodes is a list of all x positions
    for x_node in nodes:
        rho_total.append(getRefPropValue('D', 'T', T_total[-1], 'P', p_total[-1], fluid))
        try:
            rho_init = rho[-1]
            p_init = p[-1]
        except:
            rho_init = rho_total[-1]
            p_init = p_total[-1]
        node_p, node_T, node_rho, node_v, node_c, node_u = adiabatic_tube(path=path, mdot_1=mdot, p_total_2=p_total[-1],
                                                                          p_init=p_init,
                                                                          rho_total_0=rho_total[-1], rho_0=rho_init,
                                                                          T_total_0=T_total[-1], A_1=A, fluid=fluid)
        # Append the list to collect all parameters
        p.append(node_p)
        T.append(node_T)
        rho.append(node_rho)
        v.append(node_v)
        c.append(node_c)
        u.append(node_u)

        eta.append(getRefPropValue('V', 'T', T[-1], 'P', p[-1], fluid))
        Re.append(rho[-1] * v[-1] * Dh / eta[-1])
        Lambda.append(pipe_friction_factor(Re[-1], Di, Dh, k))
        dp_total.append(Lambda[-1] * (node_length / Dh) * (rho[-1] / 2) * pow(v[-1], 2))
        Ma.append(v[-1] / c[-1])

        # Values for next iteration
        p_total.append(p_total[-1] - dp_total[-1])
        T_total.append(T_total[-1])
        h_total.append(h_total[-1])

    # Remove last value from lists because some are to long
    p_total = p_total[:-1]
    T_total = T_total[:-1]
    h_total = h_total[:-1]

    # Collect all flow parameters for a DataFrame
    tube_flow = pd.DataFrame(
        {'x_nodes': nodes, 'h_total': h_total, 'u': u, 'T_total': T_total, 'T': T, 'p_total': p_total,
         'p': p, 'dp_total': dp_total, 'rho_total': rho_total, 'rho': rho, 'c': c, 'eta': eta,
         'Lambda': Lambda, 'v': v, 'Re': Re, 'Ma': Ma})

    return tube_flow

def adiabatic_tube(path, mdot_1, p_total_2, p_init, rho_total_0, rho_0, T_total_0, A_1, fluid):
    """ Calculate the parameters for the adiabatic tube element.

    :param path: path to store the debug file
    :param mdot_1: mass flow [kg/s]
    :param p_total_2: total pressure at the end of the tube [Pa]
    :param rho_total_0: total density [kg/m^3]
    :param T_total_0: total temperature [K]
    :param A_1: cross section area [m^2]
    :param fluid: fluid of the current tube
    :return: pressure, temperature, density, velocity, speed of sound, internal energy
    :type mdot_1: float
    :type p_total_2: float
    :type rho_total_0: float
    :type T_total_0: float
    :type A_1: float
    :type fluid: str
    :rtype: (float, float, float, float, float, float)

    """

    # Use total density as initial value for the local static density, same for local static pressure
    # TODO: Maybe these are bad starting values?
    # rho_iter_init = rho_total_0
    # p = p_total_2
    rho_iter_init = rho_0
    p = p_init

    # Set the total density delta to an initial value. This should converge during the iteration
    delta = RHO_ITERATION_ACCURACY * 1.1
    # Set the iteration counter to its initial values. This will prevent an endless loop.
    iter_count = 0

    # Iteration to calculate the local static density at the end of the tube (position 2)
    while abs(delta) > RHO_ITERATION_ACCURACY and iter_count < ITERATION_MAX:

        # Debug
        with open(os.path.join(path, 'debug_iteration.txt'), 'a+') as debug:
            try:
                debug.write(
                    '\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(iter_count, rho_iter_init, p, A_1 * 1e6, T, v,
                                                                          c, kappa, cp, cv))
                R = None
            except:
                debug.write('\n{}\t{}\t{}'.format(iter_count, rho_iter_init, p))

        # print('adiabatic iter = {}'.format(iter_count))

        # kappa = getRefPropValue('ISENTROPIC_EXPANSION_COEFFICIENT', 'D', rho_iter_init, 'P', p, fluid)
        cp = getRefPropValue('Cpmass', 'D', rho_iter_init, 'P', p, fluid)
        cv = getRefPropValue('Cvmass', 'D', rho_iter_init, 'P', p, fluid)
        kappa = cp / cv  # Assume an ideal gas here

        # This is just extra error checking because sometimes, RefProp does stupid things
        # In theory, this should never happen
        if kappa < 0:
            print('No - this is wrong')
            exit()

        v = mdot_1 / (rho_iter_init * A_1)

        # Try to determine the speed of sound in the current status
        # If RefProp fails because the speed of sound is not defined for a two-phase-state, assume an ideal gas, where
        # c = sqrt(kappa * R * T)
        # This could be wrong?
        # try:
        # how = 'refprop___'
        c = getRefPropValue('A', 'D', rho_iter_init, 'P', p, fluid)
        Ma = v / c
        T = T_total_0 / (1 + (kappa - 1) / 2 * pow(Ma, 2))
        # except:
        #     how = 'ideal_gas'
        #     print(
        #         'Warning - RefProp failing, assume ideal gas for {} at p = {:3.2}\u2009bar and rho = {:4.2}\u2009kg/m^3'
        #             .format(fluid, p * 1e-5, rho_iter_init))
        #     if fluid == 'Oxygen':
        #         R = GAS_CONSTANT / MOLAR_MASS_OXYGEN  # [J/(kg K)]
        #     elif fluid == 'Hydrogen' or fluid == 'H2':
        #         R = GAS_CONSTANT / MOLAR_MASS_HYDROGEN  # [J/(kg K)]
        #
        #     T = T_total_0 - ((kappa - 1) / 2) * (pow(v, 2) / (kappa * R))
        #     c = sqrt(kappa * R * T)
        #     Ma = v / c

        rho_iter_new = rho_total_0 / (1 + (kappa - 1) / 2 * pow(Ma, 2)) ** (1 / (kappa - 1))

        # Using the calculated temperature and density, a new value for the local static pressure can be determined
        p = getRefPropValue('P', 'T', T, 'D', rho_iter_new, fluid)

        # Determine the convergence, update the initial density value for the next iteration step, updating the
        # iteration counter
        delta = rho_iter_new - rho_iter_init
        rho_iter_init = rho_iter_new
        iter_count += 1

    # Print an iteration waring, if necessary
    iteration_max_warning(iter_count, 'adiabatic tube')

    # Take values from the last iteration to calculate the final parameters
    rho = rho_iter_new
    v = mdot_1 / (rho * A_1)
    c = getRefPropValue('A', 'D', rho, 'P', p, fluid)
    u = getRefPropValue('U', 'T', T, 'P', p, fluid)

    return p, T, rho, v, c, u

def pipe_friction_factor(Re, Di, Dh, k):
    """ Calculates a pipe friction coefficient.

    :param Re: Reynolds number [-]
    :param Di: inner diameter [m]
    :param Dh: hydraulic diameter [m]
    :param k: surface roughness [m]
    :return: pipe friction coefficient Lambda [-]
    :type Re: float
    :type Di: float
    :type Dh: float
    :type k: float
    :rtype: float
    """

    # Initial value and accuracy for the pipe friction factor
    LAMBDA_INIT = 0.02
    LAMBDA_ACCURACY = 1e-5

    # Determine if the geometry is an annulus or a tube
    if Di > 0:
        geometry = 'annulus'
    else:
        geometry = 'tube'

    if Re <= 2300:
        # laminar flow
        if geometry == 'annulus':
            Lambda = 95 / Re
        elif geometry == 'tube':
            Lambda = 64 / Re

    elif Re > 2300:
        # turbulent flow
        # Use the Colebrook formula to calculate the friction coefficient
        # This has to be iterative
        Lambda_init = LAMBDA_INIT
        Lambda = 1 / pow(2 * log10(2.51 / (Re * sqrt(LAMBDA_INIT)) + 0.27 / (Dh / k)), 2)

        while abs(Lambda - Lambda_init) > LAMBDA_ACCURACY:
            Lambda_init = Lambda
            Lambda = 1 / pow(2 * log10(2.51 / (Re * sqrt(Lambda_init)) + 0.27 / (Dh / k)), 2)

    return Lambda

def iteration_max_warning(iter_count, function):
    """ Prints a warning if the iteration stopped without reaching the defined accuracy.

    :param iter_count: iteration counter
    :param function: specifies the function where the warning is raised
    :return: If the iteration has converged (True) or not (False)
    :type iter_count: int
    :type function: str
    :rtype: bool
    """
    if not iter_count < ITERATION_MAX:
        print('Warning - Flow solver iteration did not reach defined accuracy in {}.'.format(function))
        # Did not reach convergence
        return False
    else:
        # Did reach convergence
        return True

def plot_pressure(file, title, x, p, p_total, p_0, p_end):
    """ Plots the local static pressure (p) and the local total pressure (p_total) over the injector length

    """
    # Assign x-value for manifold
    x = x.to_list()
    x[0] = x[1] - abs(x[1] - x[2])
    x = [node * 1e3 for node in x]

    # Plot, including conversion to bar
    plt.figure()
    plt.plot(x, [p * 1e-5 for p in p.to_list()], label='p')
    plt.plot(x, [p_t * 1e-5 for p_t in p_total.to_list()], label='p_total')
    xlim_min, xlim_max = plt.xlim()
    plt.hlines(p_0 * 1e-5, xmin=xlim_min, xmax=xlim_max, linestyles=':', color='black')
    plt.hlines(p_end * 1e-5, xmin=xlim_min, xmax=xlim_max, linestyles=':', color='black')
    plt.legend(loc='upper right')
    plt.xlabel('x [mm]')
    plt.ylabel('local pressure [bar]')
    plt.grid()
    plt.title(title)
    plt.savefig(file)

def plot_density(file, title, x, rho, rho_total):
    """ Plots the local density (rho) and the total density (rho_total) over the injector length

    """
    # Assign x-value for manifold
    x = x.to_list()
    x[0] = x[1] - abs(x[1] - x[2])
    x = [node * 1e3 for node in x]

    # Plot, including conversion to bar
    plt.figure()
    plt.plot(x, rho.to_list(), label='rho')
    plt.plot(x, rho_total.to_list(), label='rho_total')
    plt.legend(loc='lower right')
    plt.xlabel('x [mm]')
    plt.ylabel('density [kg/m^3]')
    plt.grid()
    plt.title(title)
    plt.savefig(file)

def plot_temperature(file, title, x, T, T_total, T_in):
    """ Plots the local temperature (p) and the local total temperature (T_total) over the injector length

    """
    # Assign x-value for manifold
    x = x.to_list()
    x[0] = x[1] - abs(x[1] - x[2])
    x = [node * 1e3 for node in x]

    # Plot, including conversion to bar
    plt.figure()
    plt.plot(x, T.to_list(), label='T')
    plt.plot(x, T_total.to_list(), label='T_total')
    xlim_min, xlim_max = plt.xlim()
    plt.hlines(T_in, xmin=xlim_min, xmax=xlim_max, linestyles=':', color='black')
    plt.legend(loc='lower right')
    plt.xlabel('x [mm]')
    plt.ylabel('temperature [K]')
    plt.grid()
    plt.title(title)
    plt.savefig(file)

def plot_velocity(file, title, x, v):
    """ Plots the local velocity over the injector length

    """
    # Assign x-value for manifold
    x = x.to_list()
    x[0] = x[1] - abs(x[1] - x[2])
    x = [node * 1e3 for node in x]

    # Plot velocity
    plt.figure()
    plt.plot(x, v.to_list(), label='v')
    plt.legend(loc='lower right')
    plt.xlabel('x [mm]')
    plt.ylabel('velocity [m/s]')
    plt.grid()
    plt.title(title)
    plt.savefig(file)

def plot_speed_of_sound(file, title, x, c):
    """ Plots the local speed of sound over the injector length

    """
    # Assign x-value for manifold
    x = x.to_list()
    x[0] = x[1] - abs(x[1] - x[2])
    x = [node * 1e3 for node in x]

    # Plot velocity
    plt.figure()
    plt.plot(x, c.to_list(), label='c')
    plt.legend(loc='lower right')
    plt.xlabel('x [mm]')
    plt.ylabel('speed of sound [m/s]')
    plt.grid()
    plt.title(title)
    plt.savefig(file)

def plot_geometry(file, title, x, Da, Di):
    """ Plots the geometry of the current injector.

    """
    r_a = [r * 1e3 / 2 for r in Da.to_list()]  # Outer radius [mm]
    r_i = [r * 1e3 / 2 for r in Di.to_list()]  # Inner radius [mm]

    plt.figure()
    plt.fill_between(x.to_list(), r_a, color='blue')
    plt.fill_between(x.to_list(), [r * -1 for r in r_a], color='blue')
    plt.fill_between(x.to_list(), r_i, color='white')
    plt.fill_between(x.to_list(), [r * -1 for r in r_i], color='white')
    plt.grid()
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')
    plt.title(title)
    plt.savefig(file)

# if __name__ == "__main__":
#     flow = {'h_total': [0], 'T_total': [0]}
#     flow = pd.DataFrame(flow)
#     flow = flow.append({'h_total': 12, 'T_total': 123}, ignore_index=True)
#     c = flow
#     k = flow.append(c, ignore_index=True)

def elements_to_file(path, name, elements):
    """
    :type elements: list[FlowElement]
    :type name: str
    :type path: os.path
    :rtype: None
    """

    with open(os.path.join(path, '{}_input.txt'.format(name)), 'w+') as file:
        file.write('{:<15}\tDa\tDi\tDh\tA\tL\tV\tk\tzeta\tnumber_of_elements\n'.format('type'))
        for element in elements:
            file.write('{:<15}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{}\n'.format(element.type,
                                                                       element.Da * 1e3, element.Di * 1e3, element.Dh * 1e3,
                                                                       element.A * 1e6, element.L * 1e3, element.V * 1e9, element.k * 1e3,
                                                                       float(element.zeta), element.number_of_elements))
