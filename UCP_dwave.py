# Import necessary libraries 
from dwave.system import LeapHybridCQMSampler
import dimod
from dimod import ConstrainedQuadraticModel, Binary, Real
import json
import copy
import random
import matplotlib.pyplot as plt
import numpy as np
import csv


# Load data
data_file = "rts_gmlc/2020-11-25.json"
print('loading data')
data = json.load(open(data_file, 'r'))

#INPUTS:
random.seed(19)

# How many HOURS do you want in the problem?
HOURS = 48

# How much do you want to reduce generator capasity and demand?
reduction_percentage = 0.95
# Extract data for generators and time periods
thermal_gens = data['thermal_generators']
renewable_gens = data['renewable_generators']
time_periods = {t+1 : t for t in range(HOURS)}
gen_startup_categories = {g : list(range(0, len(gen['startup']))) for (g, gen) in thermal_gens.items()}
gen_pwl_points = {g : list(range(0, len(gen['piecewise_production']))) for (g, gen) in thermal_gens.items()}

demand = data['demand']
reserves = data['reserves']

################################################################################################
def reduce_generators(thermal_gens, renewable_gens, demand, reserves, reduction_percentage, HOURS):

    # Deep copy to avoid modifying the original data
    thermal_gens = copy.deepcopy(thermal_gens)
    renewable_gens = copy.deepcopy(renewable_gens)
    demand = copy.deepcopy(demand)
    reserves = copy.deepcopy(reserves)
    
    # Calculate the number of generators to remove
    num_thermal = len(thermal_gens)
    num_renewable = len(renewable_gens)
    
    thermal_to_remove = round(num_thermal * reduction_percentage)
    renewable_to_remove = round(num_renewable * reduction_percentage)
    
    # Ensure we don't attempt to remove more generators than exist
    thermal_to_remove = min(thermal_to_remove, num_thermal)
    renewable_to_remove = min(renewable_to_remove, num_renewable)
    
    # Remove thermal generators randomly
    removed_thermal_capacity = 0
    if thermal_to_remove > 0:
        thermal_gen_keys = list(thermal_gens.keys())
        thermal_removed = random.sample(thermal_gen_keys, thermal_to_remove)
        # Sum the capacities of removed thermal generators
        removed_thermal_capacity = sum(thermal_gens[g]['power_output_maximum'] for g in thermal_removed)
        # Remove the generators from the dictionary
        for g in thermal_removed:
            del thermal_gens[g]

    # Remove renewable generators randomly
    removed_renewable_capacity = [0] * HOURS
    if renewable_to_remove > 0:
        renewable_gen_keys = list(renewable_gens.keys())
        renewable_removed = random.sample(renewable_gen_keys, renewable_to_remove)
        # Sum the per-hour capacities of removed renewable generators
        for w in renewable_removed:
            gen_capacity = renewable_gens[w]['power_output_maximum'][:HOURS]
            for t in range(HOURS):
                removed_renewable_capacity[t] += gen_capacity[t]
        # Remove the generators from the dictionary
        for w in renewable_removed:
            del renewable_gens[w]

    # Calculate total capacities before removal
    total_initial_thermal_capacity = sum(gen['power_output_maximum'] for gen in thermal_gens.values()) + removed_thermal_capacity
    total_initial_renewable_capacity = sum(sum(gen['power_output_maximum'][:HOURS]) for gen in renewable_gens.values()) + sum(removed_renewable_capacity)
    
    # Calculate total scale factor
    total_removed_capacity = removed_thermal_capacity * HOURS + sum(removed_renewable_capacity)
    total_initial_capacity = total_initial_thermal_capacity * HOURS + total_initial_renewable_capacity
    total_scale_factor = total_removed_capacity / total_initial_capacity if total_initial_capacity > 0 else 0
    total_scale_factor = min(total_scale_factor, 1.0)  # Ensure we don't scale negatively
    # Scale demand for every hour
    scaled_demand = [demand[t] * (1 - total_scale_factor) for t in range(HOURS)]

    '''

    # Scale reserves proportionally to the scaled demand
    total_initial_demand = sum(demand)
    total_scaled_demand = sum(scaled_demand)
    reserve_scale_factor = total_scaled_demand / total_initial_demand if total_initial_demand > 0 else 1
    print(f"Reserve scale factor: {reserve_scale_factor}")

    scaled_reserves = [r * reserve_scale_factor for r in reserves]
    '''

    scale_factor_reserves = removed_thermal_capacity / total_initial_thermal_capacity if total_initial_thermal_capacity > 0 else 0

    scaled_reserves = [r * (1 - scale_factor_reserves) for r in reserves]

    return thermal_gens, renewable_gens, scaled_demand, scaled_reserves

total_thermal_capacity_before = sum(gen['power_output_maximum'] for gen in thermal_gens.values())*HOURS
print(f"Total thermal capacity before function: {total_thermal_capacity_before} MW")
total_renewable_capacity_before = sum(sum(gen['power_output_maximum'][:HOURS]) for gen in renewable_gens.values())
print(f"Total renewable capacity before function: {total_renewable_capacity_before} MW")
total_demand_before = sum(demand[t] for t in range(HOURS))
print(f"Total demand before function: {total_demand_before} MW")
total_reserves_before = sum(reserves)
print(f"Total reserves before function: {total_reserves_before} MW")

thermal_gens, renewable_gens, demand, reserves = reduce_generators(thermal_gens, renewable_gens, demand, reserves, reduction_percentage, HOURS)


total_thermal_capacity_after = sum(gen['power_output_maximum'] for gen in thermal_gens.values())*HOURS
print(f"Total thermal capacity after function: {total_thermal_capacity_after} MW")
total_renewable_capacity_after = sum(sum(gen['power_output_maximum'][:HOURS]) for gen in renewable_gens.values())
print(f"Total renewable capacity after function: {total_renewable_capacity_after} MW")
total_demand_after = sum(demand)
print(f"Total demand after function: {total_demand_after} MW")
total_reserves_after = sum(reserves)
print(f"Total reserves after function: {total_reserves_after} MW")

num_thermal_gens = len(thermal_gens)
num_renewable_gens = len(renewable_gens)
print(f"Number of thermal generators: {num_thermal_gens}")
print(f"Number of renewable generators: {num_renewable_gens}")

################################################################################################


print('Building model...')
cqm = ConstrainedQuadraticModel()

cg = {}  # Generation cost
pg = {}  # Power generation above minimum
rg = {}  # Reserves
pw = {}  # Power generation
ug = {}  # Unit on/off
vg = {}  # Unit startup
wg = {}  # Unit shutdown
dg = {}  
lg = {}

for g, gen in thermal_gens.items():
    for t in time_periods:
        # Continuous variables
        cg[g, t] = dimod.Real(f'cg_{g}_{t}')
        pg[g, t] = dimod.Real(f'pg_{g}_{t}', lower_bound=0)
        rg[g, t] = dimod.Real(f'rg_{g}_{t}', lower_bound=0)
        # Binary variables
        ug[g, t] = dimod.Binary(f'ug_{g}_{t}')
        vg[g, t] = dimod.Binary(f'vg_{g}_{t}')
        wg[g, t] = dimod.Binary(f'wg_{g}_{t}')

        # Piecewise linear variables
        for l in gen_pwl_points[g]:
            lg[g, l, t] = dimod.Real(
                f'lg_{g}_{l}_{t}', lower_bound=0, upper_bound=1)

        # Startup categories
        for s in gen_startup_categories[g]:
            dg[g, s, t] = dimod.Binary(f'dg_{g}_{s}_{t}')

for w, gen in renewable_gens.items():
    for t, t_idx in time_periods.items():
        pw[w, t] = dimod.Real(('pw',w,t),lower_bound=gen['power_output_minimum'][t_idx],upper_bound=gen['power_output_maximum'][t_idx]) #This includes constraint (24)

# Objective Function
cqm.set_objective(sum(
            sum(cg[g,t] + gen['piecewise_production'][0]['cost']*ug[g,t]
                + sum( gen_startup['cost']*dg[g,s,t] for (s, gen_startup) in enumerate(gen['startup']))
            for t in time_periods)
        for g, gen in thermal_gens.items() )
        )


# Constraints
for t, t_idx in time_periods.items():
    cqm.add_constraint(sum(pg[g, t] + thermal_gens[g]['power_output_minimum'] * ug[g, t] for g in thermal_gens) + sum(pw[w, t] for w in renewable_gens) == demand[t_idx], label=f'demand_constraint_{t}') #(2)
    cqm.add_constraint(sum(rg[g, t] for g in thermal_gens) >= reserves[t_idx], label=f'reserve_constraint_{t}') #(3)
    

for g, gen in thermal_gens.items():
    if gen['unit_on_t0'] == 1:
        if gen['time_up_minimum'] - gen['time_up_t0'] >= 1:
            cqm.add_constraint(sum((ug[g, t] - 1) for t in range(1, min(gen['time_up_minimum'] - gen['time_up_t0'], HOURS)+1)) == 0, label=f'uptime_{g}_t0') #(4)
    elif gen['unit_on_t0'] == 0:
        if gen['time_down_minimum'] - gen['time_down_t0'] >= 1:
            cqm.add_constraint(sum(ug[g,t] for t in range(1, min(gen['time_down_minimum'] - gen['time_down_t0'], HOURS)+1)) == 0, label=f'downtime_{g}_t0') #(5)
    else:
        raise Exception('Invalid unit_on_t0 for generator {}, unit_on_t0={}'.format(g, gen['unit_on_t0']))

    cqm.add_constraint(ug[g,1] - gen['unit_on_t0'] -vg[g, 1] +wg[g, 1] == 0, label=f'logical_{g}_t0') #(6)

    startup_expr = sum( 
                        sum( dg[g,s,t] 
                                for t in range(
                                                max(1, gen['startup'][s+1]['lag']-gen['time_down_t0']+1),
                                                min(gen['startup'][s+1]['lag']-1, HOURS)+1
                                              )
                            )
                       for s,_ in enumerate(gen['startup'][:-1])) ## all but last
    if isinstance(startup_expr, int):
        pass
    else:
        cqm.add_constraint(startup_expr == 0, label=f'startup_allowed_{g}_t0') #(7)

    cqm.add_constraint(pg[g, 1] + rg[g, 1] - gen['unit_on_t0']*(gen['power_output_t0'] - gen['power_output_minimum']) <= gen['ramp_up_limit'], label=f'ramp_up_{g}_t0') #(8)

    cqm.add_constraint(gen['unit_on_t0']*(gen['power_output_t0'] - gen['power_output_minimum'])- pg[g, 1] <= gen['ramp_down_limit'], label=f'ramp_down_{g}_t0') #(9)

    shutdown_constr = gen['unit_on_t0']*(gen['power_output_t0']-gen['power_output_minimum']) <= gen['unit_on_t0']*(gen['power_output_maximum'] - gen['power_output_minimum']) - max((gen['power_output_maximum'] - gen['ramp_shutdown_limit']),0)*wg[g,1] #(10)

    if isinstance(shutdown_constr, bool):
        pass
    else:
        cqm.add_constraint(shutdown_constr, label=f'shutdown_{g}_t0')

for g, gen in thermal_gens.items():
    for t in time_periods:
        cqm.add_constraint(ug[g,t] >= gen['must_run'], label=f'must_run_{g}_{t}') #(11)

        if t > 1:
            cqm.add_constraint(ug[g,t] - ug[g,t-1] - (vg[g,t] - wg[g,t]) == 0, label=f'logical_{g}_{t}') #(12)

        UT = min(gen['time_up_minimum'], HOURS)
        if t >= UT:
            cqm.add_constraint(sum(vg[g,t] for t in range(t-UT+1, t+1)) - ug[g,t]<= 0, label=f'uptime_{g}_{t}') #(13)
        DT = min(gen['time_down_minimum'], HOURS)
        if t >= DT:
            cqm.add_constraint(sum(wg[g,t] for t in range(t-DT+1, t+1)) - (1-ug[g,t]) <= 0, label=f'downtime_{g}_{t}') #(14)
        cqm.add_constraint(vg[g,t] - sum(dg[g,s,t] for s,_ in enumerate(gen['startup'])) == 0, label=f'startup_select_{g}_{t}') #(16)

        cqm.add_constraint(pg[g,t]+rg[g,t] - ((gen['power_output_maximum'] - gen['power_output_minimum'])*ug[g,t] - max((gen['power_output_maximum'] - gen['ramp_startup_limit']),0)*vg[g,t]) <= 0, label=f'gen_limit1_{g}_{t}') #(17)

        if t < len(time_periods): 
            cqm.add_constraint(pg[g,t]+rg[g,t] - ((gen['power_output_maximum'] - gen['power_output_minimum'])*ug[g,t] - max((gen['power_output_maximum'] - gen['ramp_shutdown_limit']),0)*wg[g,t+1]) <= 0, label=f'gen_limit2_{g}_{t}') #(18)

        if t > 1:
            cqm.add_constraint(pg[g,t]+rg[g,t] - pg[g,t-1] <= gen['ramp_up_limit'], label=f'ramp_up_{g}_{t}') #(19)
            cqm.add_constraint(pg[g,t-1] - pg[g,t] <= gen['ramp_down_limit'], label=f'ramp_down_{g}_{t}') #(20

        piece_mw1 = gen['piecewise_production'][0]['mw']
        piece_cost1 = gen['piecewise_production'][0]['cost']
        cqm.add_constraint(pg[g,t] - sum((piece['mw'] - piece_mw1)*lg[g,l,t] for l,piece in enumerate(gen['piecewise_production'])) == 0, label=f'power_select_{g}_{t}') #(21)
        cqm.add_constraint(cg[g,t] - sum((piece['cost'] - piece_cost1)*lg[g,l,t] for l,piece in enumerate(gen['piecewise_production'])) == 0, label=f'cost_select_{g}_{t}') #(22)
        cqm.add_constraint(ug[g,t] - sum(lg[g,l,t] for l,_ in enumerate(gen['piecewise_production'])) == 0, label=f'on_select_{g}_{t}') #(23)

for g, gen in thermal_gens.items():
    for s,_ in enumerate(gen['startup'][:-1]): ## all but last
        for t in time_periods:
            if t >= gen['startup'][s+1]['lag']:
                cqm.add_constraint(dg[g,s,t] - sum(wg[g,t-i] for i in range(gen['startup'][s]['lag'], gen['startup'][s+1]['lag'])) <= 0, label=f'startup_allowed_{g}_{s}_{t}') #(15)

print("Model setup complete.")
num_variables = len(cqm.variables)
print(f"Number of variables: {num_variables}")
num_constraints = len(cqm.constraints)
print(f"Number of constraints: {num_constraints}")

# Solve the CQM
print("Solving...")

TIME_LIMIT = 150
sampler = LeapHybridCQMSampler()
#solution = sampler.sample_cqm(cqm, label="unit commitment")
solution = sampler.sample_cqm(cqm, time_limit=TIME_LIMIT, label="Unit commitment")
print("Finished sampling...")

# Total number of solutions
total_solutions = len(solution)
print(f"Total number of solutions found: {total_solutions}")

# Filter for feasible solutions
feasible_sampleset = solution.filter(lambda d: d.is_feasible)
num_feasible_solutions = len(feasible_sampleset)
print(f"Number of feasible solutions: {num_feasible_solutions}")

objective_function_values = list(solution.record["energy"])
feasible_solutions = list(feasible_sampleset.record["energy"])

# Check if any feasible solution is available
if num_feasible_solutions > 0:
    best_solution = feasible_sampleset.first
    best_solution_energy = best_solution.energy
    print("Objective function value (best feasible):", best_solution_energy)
else:
    print("No feasible solution found.")
    best_solution = solution.first
    best_solution_energy = best_solution.energy
    print("Best objective function value (from infeasible solutions):", best_solution_energy)


'''
import csv
file_path="PLOTS.csv"

new_row = [total_solutions, num_feasible_solutions, objective_function_values, feasible_solutions, best_solution_energy, TIME_LIMIT, HOURS, reduction_percentage, num_thermal_gens, num_renewable_gens, num_variables, num_constraints]

with open(file_path, mode='a', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(new_row)

print("Row added successfully!")
'''