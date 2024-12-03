from pyomo.environ import *

def generate_benders_cut(master, dual_values, theta_j, master_var_values):

    LHS = master.beta
    RHS = theta_j

    
    for constraint, dual_dict in dual_values.items(): 
        for (g, t), dual_value in dual_dict.items():
            if abs(dual_value) < 1e-6:
                continue
            if constraint == "fixing_ug":
                RHS += dual_value * (master.ug[g, t] - master_var_values["ug"][g, t])
            elif constraint == "fixing_vg":
                RHS += dual_value * (master.vg[g, t] - master_var_values["vg"][g, t])
            elif constraint == "fixing_wg":
                RHS += dual_value * (master.wg[g, t] - master_var_values["wg"][g, t])

    # Add the Benders cut to the master problem
    master.benders_cuts.add(LHS >= RHS)



'''
for (constraint, (g, t)), dual_value in dual_values.items():
    if constraint == "fixing_ug":
        if abs(dual_value) < 1e-6:
            continue
        RHS += dual_value * (master.ug[g, t] - master_solution["ug"][g, t])
        #print(dual_value)
    elif constraint == "fixing_vg":
        if abs(dual_value) < 1e-6:
            continue
        RHS += dual_value * (master.vg[g, t] - master_solution["vg"][g, t])
        #print(dual_value)
    elif constraint == "fixing_wg":
        if abs(dual_value) < 1e-6:
            continue
        RHS += dual_value * (master.wg[g, t] - master_solution["wg"][g, t])
        print(dual_value)
'''

