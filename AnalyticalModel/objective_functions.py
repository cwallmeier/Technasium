from doublet import *


def opt_t(rate, t_breakthrough_target=30, **kwargs):
    db = geoth_doub(q_m3_h=rate, t_breakthrough_target=t_breakthrough_target, **kwargs)
    db.compute_all()
    print(db.r_h, '\n', db.t_cold_parallel_yrs, t_breakthrough_target, db.q_m3_h, db.P_doublet_MW)
    return abs(db.t_breakthrough_target - db.t_cold_parallel_yrs)


def opt_power(rate, P_doublet_target=15, **kwargs):
    db = geoth_doub(q_m3_h=rate, P_doublet_target=P_doublet_target, **kwargs)
    # db.r_w = 20
    db.compute_all()
    print('const_power func, power value:%s' % P_doublet_target)
    print('\t', db.r_h, db.t_cold_parallel_yrs, db.q_m3_h, db.P_doublet_MW, P_doublet_target)
    return abs(db.P_doublet_target - db.P_doublet_MW)

# def const_power(rate, P_doublet_target=13, **kwargs):
#     db = geoth_doub(q_m3_h=rate, P_doublet_target=P_doublet_target, **kwargs)
#     db.compute_all()
#     print('const_power func, power value:%s'%P_doublet_target)
#     print('\t',db.r_h, db.t_cold_yrs, db.q_m3_h, db.P_doublet_MW, P_doublet_target)
#     return (abs(db.P_doublet_target - db.P_doublet_MW))
#
def const_t(rate,t_breakthrough_target=30, **kwargs):
    db = geoth_doub(q_m3_h=rate, **kwargs)
    db.compute_all()
    print('const_t func, t_breakthrough_target value:%s' % t_breakthrough_target)
    return(abs(db.t_breakthrough_target - db.t_cold_parallel_yrs))
