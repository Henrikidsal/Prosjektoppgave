from pyomo.environ import *

def generate_benders_cut(master, dual_values, sub_cost, master_solution):
    
    # Initialize LHS and RHS
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


'''
def generate_benders_cut(master, thermal_gens, time_periods, dual_values, sub_cost, master_solution):


    lhs_expr = master.theta
    rhs_expr = sub_cost

    for g in thermal_gens.keys():
        for t in time_periods.keys():
            constraint_name = 'fixing_ug'
            index = (g, t)
            dual = dual_values.get((constraint_name, index), 0)
            if abs(dual) < 1e-6:
                continue

            lhs_expr -= dual * master.ug[g, t]
            rhs_expr -= dual * master_solution['ug'][g, t]

    for g in thermal_gens.keys():
        for t in time_periods.keys():
            constraint_name = 'fixing_vg'
            index = (g, t)
            dual = dual_values.get((constraint_name, index), 0)
            if abs(dual) < 1e-6:
                continue

            lhs_expr -= dual * master.vg[g, t]
            rhs_expr -= dual *  master_solution['vg'][g, t]

    for g in thermal_gens.keys():
        for t in time_periods.keys():
            constraint_name = 'fixing_wg'
            index = (g, t)
            dual = dual_values.get((constraint_name, index), 0)
            if abs(dual) < 1e-6:
                continue

            lhs_expr -= dual * master.wg[g, t]
            rhs_expr -= dual * master_solution['wg'][g, t]

    # Add the Benders optimality cut to the master problem
    master.benders_cuts.add(lhs_expr >= rhs_expr)
    '''