# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:28:23 2023

@author: adaniilidis, cwallmeier
"""

import pandas as pd
import numpy as np
from CoolProp.CoolProp import PropsSI


class Doublet:

    def __init__(self,
                 # Input parameters that are likely to change:
                 # The numbers you see here are default values. Don't change them. If you want to change something,
                 # then define it in the main file (Notebook_2026) when calling Doublet.

                 # Geometry
                 depth_m=2283,  # top reservoir depth in meters
                 thickn_m=100,  # reservoir thickness in meters

                 # Geology
                 poro_percent=25,  # reservoir posority in percent
                 perm_mD=300,  # reservoir permeability mD
                 Temp_grad=30,  # How fast the temperature increases with depth (°C/km)

                 # Engineering
                 flow_m3_h=200,  # volume flow rate m3/h
                 w_space=1200,  # well spacing in meters
                 T_inj=30,  # injection temperature degC

                 # Efficiency
                 pump_effi=60,  # efficiency of the pump, in percent

                 # Input parameters that usually don't change so much. But feel free to vary them!
                 T_surf=10,  # Average temperature on the surface of the earth (°C)
                 rho_rock=2300,  # rock density kg/m3
                 Cp_rock=1,  # specific heat capacity of the rock kJ/(kgK)
                 rho_fluid=1080,  # Density of the water in kg/m3. Depends on the salt content, which can be very high.
                 Cp_fluid=4.180,  # specific heat capacity of the water in kJ/(kgK)
                 ):

        # Don't touch this block.
        # Input parameters Geometry
        self.r_d = depth_m
        self.r_h = thickn_m
        # Input parameters rock properties
        self.poro = poro_percent / 100
        self.rho_rock = rho_rock
        self.Cp_rock = Cp_rock
        self.perm_mD = perm_mD
        self.perm_m2 = self.perm_mD * 9.8692326671601e-13 * 1e-3  # reservoir permeability in m2
        # Input parameters fluid properties
        self.rho_fluid = rho_fluid
        self.Cp_fluid = Cp_fluid
        # Input parameters pT
        self.T_ref = T_surf
        self.T_grad = Temp_grad
        self.p_ref = 101325 # pressure at the surface of the earth in Pascal
        self.p_grad = 9792100 # pressure gradient with depth, in Pa/m
        # Input parameters wells
        self.q_m3_h = flow_m3_h
        self.q_m3_s = flow_m3_h / 3600  # seconds
        self.w_space = w_space
        self.w_diam = 0.2032
        self.T_inj = T_inj
        self.p_prod = self.p_ref + ((self.r_d + self.r_h / 2) * 1e-3 * self.p_grad)  # undisturbed p at reservoir depth
        self.T_prod = self.T_ref + ((self.r_d + self.r_h / 2) * 1e-3 * self.T_grad)  # production temperature degC
        self.mu_0 = PropsSI('viscosity', 'T', self.T_prod + 273.15, 'P', self.p_prod, 'Water')
        self.mu_inj = PropsSI('viscosity', 'T', self.T_inj + 273.15, 'P', self.p_prod, 'Water')
        # Input parameters economic
        self.pump_eta = pump_effi / 100


    # This part describes the physics that happen in a geothermal doublet. Everything here are definitions.
    # That means, here we tell the computer HOW to calculate something, for example how to calculate lifetime.
    # But it does not calculate it yet, it just knows how to do it.
    # The calculation happens when we call the function in the main file.

    def lmbda(self):
        """
        This function calculates how much faster the cold water moves than the cold temperature.
        This depends on the physical properties of the water and the rock.
        """
        mobility_lambda = self.poro * self.rho_fluid * self.Cp_fluid / (
                (1 - self.poro) * self.rho_rock * self.Cp_rock + self.rho_fluid * self.Cp_fluid * self.poro)
        return mobility_lambda

    def mu(self, r):
        """
        This function calculates the viscosity of the water (how "honey-like" it is). Viscosity changes
        with temperature, which is why it is different between the two wells.
        """
        if r < 0:
            mu = self.mu_inj
        elif r > self.w_space:
            mu = self.mu_0
        else:
            ratio = r / self.w_space
            mu = self.mu_inj * (1 - ratio) + self.mu_0 * ratio
        return mu

    def dp_wells(self):
        """
        This function calculates the pressure difference between the two wells. It depends on many factors. For example,
        when you want a high flowrate of water, you have to pump harder.
        """
        c_inj = self.mu_inj / self.perm_m2 * self.q_m3_s / (2 * np.pi * self.r_h)
        c_prd = self.mu_0 / self.perm_m2 * self.q_m3_s / (2 * np.pi * self.r_h)

        self.dp_MPa = np.log((self.w_space - self.w_diam / 2) / self.w_diam / 2) * (c_inj + c_prd) * 1e-6

        return self.dp_MPa

    def t_breakthrough(self):
        """
        This function calculates the moment when colder water starts to arrive at the production well. It depends on
        several factors. For example, when the wells are further apart, it takes longer.
        """
        self.t_cold_yrs = (self.poro / self.lmbda() *
                           (2 * np.pi * self.r_h) / self.q_m3_h * self.w_space ** 2 / 6) / (365 * 24)
        return self.t_cold_yrs

    def p_pumps(self):
        """
        This function calculates how much power the pumps consume.
        """
        self.p_pumps_MW = self.dp_wells() * self.q_m3_s / self.pump_eta
        return self.p_pumps_MW

    def p_doublet(self):
        """
        This function calculates how much power the doublet delivers.
        """
        self.P_doublet_kW = self.q_m3_s * self.rho_fluid * self.Cp_fluid * (self.T_prod - self.T_inj)
        self.P_doublet_MW = self.P_doublet_kW * 1e-3

        return self.P_doublet_MW
