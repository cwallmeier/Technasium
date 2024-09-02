import pandas as pd
from inspect import signature
from objective_functions import *
import time
import os
from pathlib import Path
from scipy.optimize import minimize

input_from_file = False
run_default = not input_from_file

if input_from_file:
    # Input is an excel file (or anything that can be read into a dataframe) with one row per case to run. Input may contain
    # values to be used in step 1 (geoth_doub) and step 2 (lcoh_compute)

    inputpath = Path("inputs/example_input_1000cases.xlsx")
    inputs = pd.read_excel(inputpath)

    t_start = time.time()
    # Split input dataframe into values to be used in step 1 (geoth_doub) and step 2 (lcoh_compute)
    params_doub = list(geoth_doub().__dict__.keys())
    params_econ = list(signature(lcoh_compute).parameters)
    inputs_doub = pd.DataFrame()
    inputs_econ = pd.DataFrame()
    for c in inputs.columns:
        if c in params_doub:
            inputs_doub[c] = inputs[c]
        elif c in params_econ:
            inputs_econ[c] = inputs[c]
        else:
            print("\nInput file contains column \n\"%s\"\nI don't know what to do with it, so I'll ignore it." % (str(c)))


    # apply the doublet calculation to each row
    def doublet_calculations(row):
        db_temp = geoth_doub(**row)
        db_temp.compute_all()
        return db_temp.CapEx, db_temp.p_pumps_MW, db_temp.P_doublet_MW


    results_doub = inputs_doub.apply(doublet_calculations, axis=1, result_type='expand')
    results_doub.columns = (['CapEx', 'p_pumps_MW', 'P_doublet_MW'])

    # and apply lcoh calculation
    inputs_econ = pd.concat([inputs_econ, results_doub], axis=1)
    lcoh = inputs_econ.apply(lambda row: lcoh_compute(
        CapEx=row['CapEx'],
        pumps_MW=row['p_pumps_MW'],
        p_doublet_MW=row['P_doublet_MW'])
                             , axis=1)

    howlong = round(time.time() - t_start, 3)
    print('\n\nLCOH:')
    print(lcoh)
    print('\n\n Took ' + str(howlong) + ' seconds for ' + str(len(inputs)) + ' cases')

if run_default:
    db = geoth_doub()
    db.compute_all(use_radial=True)
    db.lcoh = lcoh_compute(CapEx=db.CapEx,
                           pumps_MW=db.p_pumps_MW,
                           p_doublet_MW=db.P_doublet_MW)
    attrs = vars(db)
    print('\n '.join("%s: %s" % item for item in attrs.items()))

# optimize on power with time constraint
# opt_args = {'P_doublet_target': 11.58, 'r_h': 200}
# cons = ({'type': 'ineq', 'fun': const_t, 'args': (20,)})
# result = minimize(opt_power, 100, bounds=[[1, 500]], args=(20,), method='SLSQP', constraints=cons)
# opt_args = {'P_doublet_target': 20, 't_breakthrough_target': 35, 'r_h': 200}
# cons = ({'type': 'ineq', 'fun': lambda x, **kwargs: opt_t(x, **opt_args)})
# result = minimize(lambda rate: opt_power(rate, **opt_args), 100, bounds=[[1, 500]], method='SLSQP', constraints=cons)
# print(result.x[0].round())
#
# db = geoth_doub(q_m3_h=result.x[0].round(), **opt_args)
# db.compute_all()
# attrs = vars(db)
# print('\n '.join("%s: %s" % item for item in attrs.items()))
#
#
#
# # optimize on time with power constraint
# # cons = ({'type':'ineq', 'fun':const_power, 'args':(17,), opt_args})
# # result = minimize(opt_t, 40, bounds=[[1,500]], args=(20,), method='SLSQP', constraints=cons)
# opt_args = {'P_doublet_target':20, 't_breakthrough_target':35, 'r_h':200}
# cons = ({'type':'ineq', 'fun':lambda x, **kwargs: opt_power(x, **opt_args)})
# result = minimize(lambda rate: opt_t(rate, **opt_args), 40, bounds=[[1,500]], method='SLSQP', constraints=cons)
# print(result.x[0].round())
#
# db = geoth_doub(q_m3_h=result.x[0].round(), **opt_args)
# db.compute_all()
# attrs = vars(db)
# print('\n '.join("%s: %s" % item for item in attrs.items()))
