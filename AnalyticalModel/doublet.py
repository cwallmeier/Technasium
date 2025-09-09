# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:28:23 2023

@author: adaniilidis, cwallmeier
"""

import pandas as pd
import numpy as np
from CoolProp.CoolProp import PropsSI


class make_doublet():

    def __init__(self,
                 # Input parameters Geometry
                 depth_m=2283.33334,  # top reservoir depth in meters
                 thickn_m=100,  # reservoir thickness in meters
                 width_m=1200,  # reservoir width in meters

                 # Input parameters rock properties
                 poro=0.15,  # reservoir posority
                 perm_mD=300,  # reservoir permeability mD
                 rho_rock=2300,  # rock density kg/m3
                 Cp_rock=1,  # 0.8532,  # specific heat capacity rock kJ/(kgK)
                 cond_conf=2, # thermal conductivity of confining layers, in W/(mK)

                 # Input parameters fluid properties
                 rho_fluid=1000,  # fluid density kg/m3
                 Cp_fluid=4.180,  # specific heat capacity fluid kJ/(kgK)

                 # Input parameters pT
                 T_surf=10,  # refrence temperature at surface (degC)
                 T_grad=30,  # temperature gradient (degC/km)
                 p_surf=101325,  # reference pressure at surface (Pa)
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
        self.r_d = depth_m
        self.r_h = thickn_m
        self.r_w = width_m
        self.A = self.r_h * self.r_w  # cross sectional area (m2)

        # Input parameters rock properties
        self.poro = poro
        self.rho_rock = rho_rock
        self.Cp_rock = Cp_rock
        self.cond_conf = cond_conf
        self.perm_mD = perm_mD
        self.perm_m2 = self.perm_mD * 9.8692326671601e-13 * 1e-3  # reservoir permeability in m2

        # Input parameters fluid properties
        self.rho_fluid = rho_fluid
        self.Cp_fluid = Cp_fluid

        # Input parameters pT
        self.T_ref = T_surf
        self.T_grad = T_grad
        self.p_ref = p_surf
        self.p_grad = p_grad

        # Input parameters wells
        self.q_m3_h = q_m3_h
        self.q_m3_s = q_m3_h / 3600  # seconds
        self.w_space = w_space
        self.w_diam = w_diam
        self.T_inj = T_inj
        self.p_prod = self.p_ref + ((self.r_d + self.r_h / 2) * 1e-3 * self.p_grad)  # undisturbed p at reservoir depth
        self.T_prod = self.T_ref + ((self.r_d + self.r_h / 2) * 1e-3 * self.T_grad)  # production temperature degC
        self.mu_0 = PropsSI('viscosity', 'T', self.T_prod + 273.15, 'P', self.p_prod, 'Water')
        self.mu_inj = PropsSI('viscosity', 'T', self.T_inj + 273.15, 'P', self.p_prod, 'Water')

        # Input parameters economic
        self.pump_eta = pump_eta

    def lmbda(self):
        mobility_lambda = self.poro * self.rho_fluid * self.Cp_fluid / (
                (1 - self.poro) * self.rho_rock * self.Cp_rock + self.rho_fluid * self.Cp_fluid * self.poro)
        return mobility_lambda

    def mu(self, r):
        """
        :param r: position along doublet-line in m, measured from injector towards producer
        :return: viscosity in Pas at that position
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

        c_inj = self.mu_inj / self.perm_m2 * self.q_m3_h / (2 * np.pi * self.r_h)
        c_prd = self.mu_0 / self.perm_m2 * self.q_m3_h / (2 * np.pi * self.r_h)

        self.dp_MPa = np.log((self.w_space - self.w_diam / 2) / self.w_diam / 2) * (c_inj + c_prd) * 1e-6

        return self.dp_MPa

    def t_breakthrough(self):
        """
        :return: arrival time of cold front at producer well
        """
        self.t_cold_yrs = (self.poro / self.lmbda() *
                           (2 * np.pi * self.r_h) / self.q_m3_h * self.w_space ** 2 / 6) / (365 * 24 * 3600)
        return self.t_cold_yrs

    def p_pumps(self):
        self.p_pumps_MW = self.dp_wells() * self.q_m3_s / self.pump_eta
        return self.p_pumps_MW

    def p_doublet(self):
        """ power produced by doublet in MW """
        self.P_doublet_kW = self.q_m3_s * self.rho_fluid * self.Cp_fluid * (self.T_prod - self.T_inj)
        self.P_doublet_MW = self.P_doublet_kW * 1e-3

        return self.P_doublet_MW

    # ------------------------  Maybe for animating: ----------------------------------------------
    def v_coldfront(self, r):
        # pressure solution is undefined (towards inf) at the well locations.
        # So change r to outer well diameter in those cases
        if r == 0:
            r = self.w_diam/2
        elif r == self.w_space:
            r = self.w_space - self.w_diam/2
        v_darcy = (self.q_m3_h / (2 * np.pi * self.r_h * r) +
                   self.q_m3_h / (2 * np.pi * self.r_h * (self.w_space - r)))
        v_cold = self.lmbda() / self.poro * v_darcy
        return v_cold
