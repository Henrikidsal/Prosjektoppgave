from pyomo.environ import *

def generate_benders_cut(master, dual_values, sub_cost, master_solution):
    
    LHS = master.theta
    RHS = sub_cost

    # Add terms from duals for ug, vg, wg
    for (constr_name, (g, t)), dual_value in dual_values.items():
        if constr_name == "fixing_ug":
            LHS -= dual_value * master.ug[g, t]
            RHS -= dual_value * master_solution["ug"][g, t]
        elif constr_name == "fixing_vg":
            LHS -= dual_value * master.vg[g, t]
            RHS -= dual_value * master_solution["vg"][g, t]
        elif constr_name == "fixing_wg":
            LHS -= dual_value * master.wg[g, t]
            RHS -= dual_value * master_solution["wg"][g, t]
            
    # Add the Benders cut to the master problem
    master.benders_cuts.add(LHS >= RHS)