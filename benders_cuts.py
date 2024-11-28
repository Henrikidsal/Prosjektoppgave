# benders_cuts.py

from pyomo.environ import *
import pyomo.environ as pyo

def generate_benders_cut(master, thermal_gens, time_periods, master_solution, subproblem_model, dual_values, solver_status):

    if not solver_status:
        print('Generating feasibility cut.')
        # Construct the feasibility cut

        # Initialize the cut expression
        cut_expr = 0

        # Loop over dual variables from the feasibility dual problem
        for (var_name, index), dual_value in dual_values.items():
            if dual_value is None or abs(dual_value) < 1e-6:
                continue

            # The dual variables are associated with the constraints fixing the master variables
            # These constraints are of the form m.ug[g,t] == master_solution['ug'][g,t], etc.
            if var_name == 'fix_ug':
                g, t = index
                cut_expr += dual_value * (master.ug[g, t] - master_solution['ug'][(g, t)])
            elif var_name == 'fix_vg':
                g, t = index
                cut_expr += dual_value * (master.vg[g, t] - master_solution['vg'][(g, t)])
            elif var_name == 'fix_wg':
                g, t = index
                cut_expr += dual_value * (master.wg[g, t] - master_solution['wg'][(g, t)])

        # Add the feasibility cut to the master problem
        master.benders_cuts.add(cut_expr <= 0)
        print('Feasibility cut added to the master problem.')

    else:
        print('Generating optimality cut.')

        # Calculate the objective value of the subproblem
        subproblem_obj = pyo.value(subproblem_model.obj)

        # Start constructing the RHS of the Benders cut
        rhs = subproblem_obj

        # Loop over dual variables
        for (constr_name, index), dual_value in dual_values.items():
            if dual_value is None or abs(dual_value) < 1e-6:
                continue

            if constr_name == 'demand':
                t = index
                for g in thermal_gens.keys():
                    gen = thermal_gens[g]
                    coeff = gen['power_output_minimum']
                    rhs += dual_value * coeff * (master.ug[g, t] - master_solution['ug'][(g, t)])

            elif constr_name == 'gen_limit1':
                g, t = index
                gen = thermal_gens[g]
                max_output = gen['power_output_maximum'] - gen['power_output_minimum']
                startup_limit = max(gen['power_output_maximum'] - gen['ramp_startup_limit'], 0)
                coeff_ug = -max_output
                coeff_vg = startup_limit
                rhs += dual_value * (coeff_ug * (master.ug[g, t] - master_solution['ug'][(g, t)])
                                     + coeff_vg * (master.vg[g, t] - master_solution['vg'][(g, t)]))

            elif constr_name == 'gen_limit2':
                g, t = index
                gen = thermal_gens[g]
                max_output = gen['power_output_maximum'] - gen['power_output_minimum']
                shutdown_limit = max(gen['power_output_maximum'] - gen['ramp_shutdown_limit'], 0)
                coeff_ug = -max_output
                coeff_wg = shutdown_limit
                if t + 1 in time_periods:
                    rhs += dual_value * (coeff_ug * (master.ug[g, t] - master_solution['ug'][(g, t)])
                                         + coeff_wg * (master.wg[g, t + 1] - master_solution['wg'][(g, t + 1)]))

            elif constr_name == 'on_select':
                g, t = index
                rhs += dual_value * ( - (master.ug[g, t] - master_solution['ug'][(g, t)]))

            elif constr_name == 'ramp_up':
                g, t = index
                if t > 1:
                    gen = thermal_gens[g]
                    min_output = gen['power_output_minimum']
                    rhs += dual_value * (-min_output * (master.ug[g, t - 1] - master_solution['ug'][(g, t - 1)]))

            elif constr_name == 'ramp_down':
                g, t = index
                if t > 1:
                    gen = thermal_gens[g]
                    min_output = gen['power_output_minimum']
                    rhs += dual_value * (min_output * (master.ug[g, t - 1] - master_solution['ug'][(g, t - 1)]))

            elif constr_name == 'rampupt0':
                g = index
                gen = thermal_gens[g]
                min_output = gen['power_output_minimum']
                rhs += dual_value * min_output * (master.ug[g, 1] - master_solution['ug'][(g, 1)])

            elif constr_name == 'rampdownt0':
                g = index
                gen = thermal_gens[g]
                min_output = gen['power_output_minimum']
                rhs += dual_value * (-min_output * (master.ug[g, 1] - master_solution['ug'][(g, 1)]))

            elif constr_name == 'power_select':
                g, t = index
                rhs += dual_value * ( - (master.ug[g, t] - master_solution['ug'][(g, t)]))

            # Include other constraints if necessary

        # Add the optimality cut to the master problem
        master.benders_cuts.add(master.theta >= rhs)
        print('Optimality cut added to the master problem.')

