from pyomo.environ import *

def generate_benders_cut(master, m, dual_values, theta_j, master_var_values):

    LHS = master.beta
    RHS = theta_j

    for constraint, dual_dict in dual_values.items(): 
        for (g, t), dual_value in dual_dict.items():
            if abs(dual_value) < 1e-6:
                continue
            if constraint == "fixing_ug":
                RHS += abs(dual_value) * (master.ug[g, t] - master_var_values["ug"][g, t])
            elif constraint == "fixing_vg":
                RHS += abs(dual_value) * (master.vg[g, t] - master_var_values["vg"][g, t])
            elif constraint == "fixing_wg":
                RHS += abs(dual_value) * (master.wg[g, t] - master_var_values["wg"][g, t])


    # Add the Benders cut to the master problem
    master.benders_cuts.add(LHS >= RHS)

