from objective_functions import *

# I want to get the different lifetimes for different spacings

# make a vector of well spacings
ws = np.arange(100, 2001, 100)

# make a dataframe where all the breakthrough times go
cols = ['parallel', 'radial', 'fullDarcy']
df_bt = pd.DataFrame(columns=cols, index=ws)

# make a dataframe where all the pressures go
cols = ['p_Inj', 'p_Prd', 'dp']
df_dp = pd.DataFrame(columns=cols, index=ws)

for s in ws:
    db = geoth_doub(w_space=s)
    db.compute_all(use_radial=True)

    df_bt.loc[s] = [db.t_cold_parallel_yrs, db.t_cold_radial_yrs, db.t_cold_dp_yrs]
    df_dp.loc[s] = [db.p_Inj_radial_MPa, db.p_Prd_radial_MPa, db.dp_radial_MPa]

df_bt.to_excel('outputs/BT_by_WS_analytical.xlsx')
df_dp.to_excel('outputs/P_by_WS_analytical.xlsx')
