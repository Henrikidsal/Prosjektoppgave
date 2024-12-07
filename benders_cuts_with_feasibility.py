# benders_cuts.py

from pyomo.environ import *
import pyomo.environ as pyo

def generate_benders_cut(sub_problem_status, master, thermal_gens, data, time_periods, master_solution, subproblem_model, dual_values):

    if sub_problem_status == True:
        print('Sub problem is feasable, adding optimality cut to the master problem')
        generate_optimality_cut(master=master, thermal_gens=thermal_gens , time_periods=time_periods, master_solution=master_solution, subproblem_model=subproblem_model, dual_values=dual_values)
    elif sub_problem_status == False:
        print('Sub problem is infeasable, adding feasability cut to the master problem')
        generate_feasibility_cut(master, thermal_gens, data, time_periods, dual_values)
    else:
        print("Something is wrong")

def generate_optimality_cut(master, thermal_gens, time_periods, master_solution, subproblem_model, dual_values):
    
    #Think maybe there is something wrong with setting up the optimality cuts, and also the connection with master and sub problem. why is the theta not given a number anywhere? and whats supproblem objective doing?
    subproblem_obj = pyo.value(subproblem_model.obj)

    lhs_expr = master.theta
    rhs_value = subproblem_obj

    for (constr_name, index), dual_value in dual_values.items():
        if dual_value is None or abs(dual_value) < 1e-6:
            continue

        if constr_name == 'demand':
            t = index
            for g in thermal_gens.keys():
                gen = thermal_gens[g]
                coeff = dual_value * gen['power_output_minimum']
                lhs_expr += coeff * master.ug[g, t]
                rhs_value += coeff * master_solution['ug'][(g, t)]

        elif constr_name == 'gen_limit1':
            g, t = index
            gen = thermal_gens[g]
            max_output = gen['power_output_maximum'] - gen['power_output_minimum']
            startup_limit = max(gen['power_output_maximum'] - gen['ramp_startup_limit'], 0)
            coeff_ug = dual_value * (-max_output)
            coeff_vg = dual_value * startup_limit
            lhs_expr += coeff_ug * master.ug[g, t] + coeff_vg * master.vg[g, t]
            rhs_value += coeff_ug * master_solution['ug'][(g, t)] + coeff_vg * master_solution['vg'][(g, t)]

        elif constr_name == 'gen_limit2':
            g, t = index
            gen = thermal_gens[g]
            max_output = gen['power_output_maximum'] - gen['power_output_minimum']
            shutdown_limit = max(gen['power_output_maximum'] - gen['ramp_shutdown_limit'], 0)
            coeff_ug = dual_value * (-max_output)
            coeff_wg = dual_value * shutdown_limit
            lhs_expr += coeff_ug * master.ug[g, t]
            rhs_value += coeff_ug * master_solution['ug'][(g, t)]
            if t + 1 in time_periods:
                lhs_expr += coeff_wg * master.wg[g, t + 1]
                rhs_value += coeff_wg * master_solution['wg'][(g, t + 1)]

        elif constr_name == 'on_select':
            g, t = index
            coeff = dual_value
            lhs_expr += coeff * master.ug[g, t]
            rhs_value += coeff * master_solution['ug'][(g, t)]

    master.benders_cuts.add(lhs_expr >= rhs_value)

def generate_feasibility_cut(master, thermal_gens, data, time_periods, dual_values):

    lhs_expr = 0

    for (constr_name, index), dual_value in dual_values.items():
        if dual_value is None or abs(dual_value) < 1e-6:
            continue  # Skip negligible duals

        if constr_name == 'demand':
            t = index  # t in time_periods
            dual_demand_t = dual_value
            t_idx = time_periods[t]
            b_i = data['demand'][t_idx]
            # A_i_x = sum_{g} Pmin_g * master.ug[g, t]
            A_i_x = sum(thermal_gens[g]['power_output_minimum'] * master.ug[g, t] for g in thermal_gens.keys())
            lhs_expr += dual_demand_t * (b_i - A_i_x)
        elif constr_name == 'reserves':
            t = index
            dual_reserve_t = dual_value
            t_idx = time_periods[t]
            b_i = data['reserves'][t_idx]
            # No master variables in reserves constraints
            lhs_expr += dual_reserve_t * b_i

    # Add the feasibility cut to the master problem
    master.benders_cuts.add(lhs_expr <= 0)









