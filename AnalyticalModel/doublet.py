# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:28:23 2023

@author: adaniilidis, cwallmeier
"""

import pandas as pd
import numpy as np
from CoolProp.CoolProp import PropsSI

class geoth_doub():

    def __init__(self,
                 # Input parameters Geometry
                 r_d=2283.33334,  # top reservoir depth in meters
                 r_h=100,  # reservoir thickness in meters
                 r_w=500,  # reservoir width in meters

                 # Input parameters rock properties
                 phi=0.15,  # reservoir posority
                 k_mD=300,  # reservoir permeability mD
                 rho_rock=2300,  # rock density kg/m3
                 Cp_rock=1,  # 0.8532,  # specific heat capacity rock kJ/(kgK)

                 # Input parameters fluid properties
                 rho_fluid=1000,  # fluid density kg/m3
                 Cp_fluid=4.180,  # specific heat capacity fluid kJ/(kgK)

                 # Input parameters pT
                 T_ref=10,  # refrence temperature at surface (degC)
                 T_grad=30,  # temperature gradient (degC/km)
                 p_ref=101325,  # reference pressure at surface (Pa)
                 p_grad=9792100,  # pressure gradient (Pa/km)

                 # Input parameters wells
                 q_m3_h=172.249,  # 250,  # volume flow rate m3/h
                 w_space=1200,  # 1200,  # well spacing in meters
                 w_diam=0.2032,  # well diameter in meters
                 T_inj=30,  # injection temperature degC

                 # Input parameters economic
                 pump_eta=0.5,  # 
                 ):
        
        # Input parameters Geometry
        self.r_d = r_d
        self.r_h = r_h
        self.r_w = r_w
        self.A = self.r_h * self.r_w  # cross sectional area (m2)

        # Input parameters rock properties
        self.phi = phi
        self.rho_rock = rho_rock
        self.Cp_rock = Cp_rock
        self.k_mD = k_mD
        self.k_m2 = k_mD * 9.8692326671601e-13 * 1e-3  # reservoir permeability in m2

        # Input parameters fluid properties
        self.rho_fluid = rho_fluid
        self.Cp_fluid = Cp_fluid

        # Input parameters pT
        self.T_ref = T_ref
        self.T_grad = T_grad
        self.p_ref = p_ref
        self.p_grad = p_grad

        # Input parameters wells
        self.q_m3_h = q_m3_h
        self.q_m3_s = q_m3_h / 3600  # seconds
        self.w_space = w_space
        self.w_diam = w_diam
        self.T_inj = T_inj
        self.p_prod = self.p_ref + ((self.r_d + self.r_h / 2) * 1e-3 * self.p_grad)  # undisturbed p at reservoir depth
        self.T_prod = self.T_ref + ((self.r_d + self.r_h / 2) * 1e-3 * self.T_grad)  # production temperature degC
        self.mu_prod = PropsSI('viscosity', 'T', self.T_prod + 273.15, 'P', self.p_prod, 'Water'),  # Pa/s
        self.mu_prod = self.mu_prod[0]
        # !!!!!!!!!!!! CHANGED PRESSURE TO P_prod for mu_inj, because why use surface pressure?
        # self.mu_inj = PropsSI('viscosity', 'T', self.T_inj + 273.15, 'P', self.p_ref, 'Water'),  # Pa/s
        self.mu_inj = float(PropsSI('viscosity', 'T', self.T_inj + 273.15, 'P', self.p_prod, 'Water')),  # Pa/s
        self.mu_inj = self.mu_inj[0]
        self.mu_avg = (self.mu_prod + self.mu_inj) / 2

        # Input parameters economic
        self.pump_eta = pump_eta



    def compute_all(self):

        self.dp_wells()
        self.v_coldfront()
        self.t_breakthrough()
        self.p_doublet()
        self.p_pumps()
        self.q_sodm()

    def dp_wells(self):
        """
        Option 1 for calculating pressure drop between Inj and Prd in MPa
        Assumes parallel flow from producer to injector, with cold front moving as a flat surface that covers the entire
        (rectancular) cross-section of the reservoir
        """
        self.dp_MPa = 1e-6 * self.q_m3_s * self.mu_avg * self.w_space / (self.k_m2 * self.A)

        return self.dp_MPa


    def p_pumps(self):

        self.p_pumps_MW = self.dp_wells() * self.q_m3_s / self.pump_eta
        return self.p_pumps_MW

    def mobility_lambda(self):
        mobility_lambda = self.phi * self.rho_fluid * self.Cp_fluid / (
                (1 - self.phi) * self.rho_rock * self.Cp_rock + self.rho_fluid * self.Cp_fluid * self.phi)
        return mobility_lambda

    def v_coldfront(self):
        """
        Option 1 for calculating velocity of the cold front in m/s.
        Assumes parallel flow from producer to injector, with cold front moving as a flat surface that covers the entire
        (rectancular) cross-section of the reservoir
        """
        self.v_darcy = self.q_m3_s / self.A  # in m/s
        self.v_cfront = self.mobility_lambda() / self.phi * self.v_darcy
        return self.v_cfront


    def t_breakthrough(self):
        """
        Option 1 for time of thermal breakthrough in years.
        Assumes parallel flow from producer to injector, with cold front moving as a flat surface that covers the entire
        (rectancular) cross-section of the reservoir
        """
        self.t_cold_yrs = self.w_space / (self.v_coldfront() * 365 * 24 * 60 * 60)

        return self.t_cold_yrs


    def p_doublet(self):
        """ power produced by doublet in MW """
        self.P_doublet_kW = self.q_m3_s * self.rho_fluid * self.Cp_fluid * (self.T_prod - self.T_inj)
        self.P_doublet_MW = self.P_doublet_kW * 1e-3

        return self.P_doublet_MW
