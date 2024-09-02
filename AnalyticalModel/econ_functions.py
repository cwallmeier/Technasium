import numpy as np
import pandas as pd


def drill_costs(drill_depth=1, well_number=2):
    drill_cost = 1480000 + 1150 * drill_depth + 0.3 * drill_depth ** 2
    return drill_cost * well_number


def lcoh_compute(CapEx=1e5,
                 pumps_MW=1,
                 p_doublet_MW=16.4,
                 elec_price_MWh=120,
                 annual_discount_rate=0.05,
                 pump_cost=5e5,
                 pump_replace=5,
                 annual_OpEx_percent_CapEx=0.05,
                 lcoh_years=30,
                 he_efficiency=0.85):
    temp_econ = pd.DataFrame(index=range(lcoh_years), columns=['CapEx'], dtype=float)
    temp_econ['timeperiods'] = range(len(temp_econ))
    temp_econ['CapEx'] = 0.
    temp_econ.loc[0, 'CapEx'] = CapEx + pump_cost
    temp_econ.loc[::pump_replace, 'CapEx'] += pump_cost
    temp_econ['OpEx'] = temp_econ['CapEx'].cumsum() * annual_OpEx_percent_CapEx
    year_hours = 24 * 365.25
    temp_econ['OpEx_pump'] = elec_price_MWh * pumps_MW * year_hours

    temp_econ['Heat produced'] = p_doublet_MW * year_hours * he_efficiency
    temp_econ['Costs total'] = temp_econ['CapEx'] + temp_econ['OpEx'] + temp_econ['OpEx_pump']
    temp_econ['Discounted costs'] = (
            temp_econ['Costs total'] / ((1 + annual_discount_rate) ** temp_econ['timeperiods'])).cumsum()
    temp_econ['Discounted energy'] = (
            temp_econ['Heat produced'] / ((1 + annual_discount_rate) ** temp_econ['timeperiods'])).cumsum()
    temp_econ['LCOH (\u20ac/MWh)'] = temp_econ['Discounted costs'] / temp_econ['Discounted energy']

    return temp_econ['LCOH (\u20ac/MWh)'].iloc[-1]
