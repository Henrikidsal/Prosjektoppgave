from pyomo.environ import *

def generate_benders_cut(master, thermal_gens, time_periods, dual_values, sub_cost, master_solution):


    lhs_expr = master.theta
    rhs_expr = sub_cost

    # Demand constraints
    for t in time_periods.keys():
        constraint_name = 'demand'
        index = t
        dual = dual_values.get((constraint_name, index), 0)
        if abs(dual) < 1e-6:
            continue

        for g, gen in thermal_gens.items():
            lhs_expr -= dual * gen['power_output_minimum'] * master.ug[g, t]
            rhs_expr -= dual * gen['power_output_minimum'] * master_solution['ug'][g, t]

    # On-Select constraints
    for g in thermal_gens.keys():
        for t in time_periods.keys():
            constraint_name = 'on_select'
            index = (g, t)
            dual = dual_values.get((constraint_name, index), 0)
            if abs(dual) < 1e-6:
                continue

            lhs_expr -= dual * master.ug[g, t]
            ug_fixed = master_solution['ug'][g, t]
            rhs_expr -= dual * ug_fixed

    # Generator Output Limits After Startup (gen_limit1)
    for g in thermal_gens.keys():
        gen = thermal_gens[g]
        for t in time_periods.keys():
            constraint_name = 'gen_limit1'
            index = (g, t)
            dual = dual_values.get((constraint_name, index), 0)
            if abs(dual) < 1e-6:
                continue

            c1 = gen['power_output_maximum'] - gen['power_output_minimum']
            c2 = max(0, gen['power_output_maximum'] - gen['ramp_startup_limit'])

            lhs_expr -= dual * (c1 * master.ug[g, t] - c2 * master.vg[g, t])

            ug_fixed = master_solution['ug'][g, t]
            vg_fixed = master_solution['vg'][g, t]

            rhs_expr -= dual * (c1 * ug_fixed - c2 * vg_fixed)

    # Generator Output Limits Before Shutdown (gen_limit2)
    for g in thermal_gens.keys():
        gen = thermal_gens[g]
        for t in time_periods.keys():
            if t == max(time_periods.keys()):  # Skip the last time period
                continue
            constraint_name = 'gen_limit2'
            index = (g, t)
            dual = dual_values.get((constraint_name, index), 0)
            if abs(dual) < 1e-6:
                continue

            c1 = gen['power_output_maximum'] - gen['power_output_minimum']
            c2 = max(0, gen['power_output_maximum'] - gen['ramp_shutdown_limit'])

            lhs_expr -= dual * (c1 * master.ug[g, t] - c2 * master.wg[g, t + 1])
            ug_fixed = master_solution['ug'][g, t]
            wg_fixed = master_solution['wg'][g, t + 1]
            rhs_expr -= dual * (c1 * ug_fixed - c2 * wg_fixed)

    # Add the Benders optimality cut to the master problem
    master.benders_cuts.add(lhs_expr >= rhs_expr)