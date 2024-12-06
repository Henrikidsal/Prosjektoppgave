# Import necessary libraries 
from dwave.system import LeapHybridCQMSampler
from dimod import ConstrainedQuadraticModel, Binary, Real
import json
import copy
import random

# Load data
data_file = "rts_gmlc/2020-02-09.json"
print('Loading data...')
data = json.load(open(data_file, 'r'))

#INPUTS:
random.seed(19)

# How many HOURS do you want in the problem?
HOURS = 24

# How much do you want to reduce generator capasity and demand?
reduction_percentage = 0.3

# Extract data for generators and time periods
data["time_periods"] = HOURS
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
    
    # Step 1: Update demand by subtracting removed renewable capacities per hour
    demand_updated_1 = []
    for t in range(HOURS):
        updated_demand = demand[t] - removed_renewable_capacity[t]
        demand_updated_1.append(max(updated_demand, 0))

    # Calculate total thermal capacity before removal
    total_initial_thermal_capacity = sum(gen['power_output_maximum'] for gen in thermal_gens.values()) + removed_thermal_capacity
    
    
    # Calculate scale_factor based on thermal capacity removed
    scale_factor = removed_thermal_capacity / total_initial_thermal_capacity

    total_scale_factor = (removed_thermal_capacity + sum(removed_renewable_capacity)) / (total_initial_thermal_capacity + total_renewable_capacity_before)
    
    # Step 2: Update demand by applying the scale_factor
    demand_updated_2 = [d * (1 - scale_factor) for d in demand_updated_1]
    
    # Update reserves by applying the scale_factor
    reserves_updated = [r * (1 - scale_factor) for r in reserves]

    demandd = [d * total_scale_factor for d in demand]

    
    return thermal_gens, renewable_gens, demandd, reserves_updated


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
ug = {}  # Unit on/off
vg = {}  # Unit startup
wg = {}  # Unit shutdown
dg = {}  
lg = {}
pw = {}  # Power generation

#This section is only here in Dwave, not neccesary in pyomo
# Define variables for thermal generators
for g, gen in thermal_gens.items():
    for t in time_periods:
        # Continuous variables
        cg[g, t] = Real(f'cg_{g}_{t}')
        pg[g, t] = Real(f'pg_{g}_{t}', lower_bound=0)
        rg[g, t] = Real(f'rg_{g}_{t}', lower_bound=0)
        # Binary variables
        ug[g, t] = Binary(f'ug_{g}_{t}')
        vg[g, t] = Binary(f'vg_{g}_{t}')
        wg[g, t] = Binary(f'wg_{g}_{t}')

        # Piecewise linear variables
        for l in gen_pwl_points[g]:
            lg[g, l, t] = Real(
                f'lg_{g}_{l}_{t}', lower_bound=0, upper_bound=1)

        # Startup categories
        for s in gen_startup_categories[g]:
            dg[g, s, t] = Binary(f'dg_{g}_{s}_{t}')

# Define variables for renewable generators
for w, gen in renewable_gens.items():
    for t in time_periods:
        pw[w, t] = Real(f'pw_{w}_{t}', lower_bound=0)

# Objective Function
objective = 0
for g, gen in thermal_gens.items():
    piecewise_cost0 = gen['piecewise_production'][0]['cost']
    for t in time_periods:
        # Generation cost and startup costs
        gen_cost = cg[g, t]
        startup_cost = sum(
            gen['startup'][s]['cost'] * dg[g, s, t]
            for s in gen_startup_categories[g]
        )
        objective += gen_cost + piecewise_cost0 * ug[g, t] + startup_cost

# Set the objective in the CQM
cqm.set_objective(objective)

# Constraints
for t, t_idx in time_periods.items():
    cqm.add_constraint(sum(pg[g, t] + thermal_gens[g]['power_output_minimum'] * ug[g, t] for g in thermal_gens) + sum(pw[w, t] for w in renewable_gens) == data['demand'][t_idx], label=f'demand_constraint_{t}') #(2)
    cqm.add_constraint(sum(rg[g, t] for g in thermal_gens) >= data['reserves'][t_idx], label=f'reserve_constraint_{t}') #(3)
    

for g, gen in thermal_gens.items():
    if gen['unit_on_t0'] == 1:
        if gen['time_up_minimum'] - gen['time_up_t0'] >= 1:
            cqm.add_constraint(sum((ug[g, t] - 1) for t in range(1, min(gen['time_up_minimum'] - gen['time_up_t0'], data['time_periods'])+1)) == 0, label=f'uptime_{g}_t0') #(4)
    elif gen['unit_on_t0'] == 0:
        if gen['time_down_minimum'] - gen['time_down_t0'] >= 1:
            cqm.add_constraint(sum(ug[g,t] for t in range(1, min(gen['time_down_minimum'] - gen['time_down_t0'], data['time_periods'])+1)) == 0, label=f'downtime_{g}_t0') #(5)
    else:
        raise Exception('Invalid unit_on_t0 for generator {}, unit_on_t0={}'.format(g, gen['unit_on_t0']))

    cqm.add_constraint(ug[g,1] - gen['unit_on_t0'] -vg[g, 1] +wg[g, 1] == 0, label=f'logical_{g}_t0') #(6)

    startup_expr = sum( 
                        sum( dg[g,s,t] 
                                for t in range(
                                                max(1, gen['startup'][s+1]['lag']-gen['time_down_t0']+1),
                                                min(gen['startup'][s+1]['lag']-1, data['time_periods'])+1
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

        UT = min(gen['time_up_minimum'], data['time_periods'])
        if t >= UT:
            cqm.add_constraint(sum(vg[g,t] for t in range(t-UT+1, t+1)) - ug[g,t]<= 0, label=f'uptime_{g}_{t}') #(13)
        DT = min(gen['time_down_minimum'], data['time_periods'])
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

for w, gen in renewable_gens.items():
    for t, t_idx in time_periods.items():
        pw[w, t].lower_bound = gen['power_output_minimum'][t_idx] #(24)
        pw[w, t].upper_bound = gen['power_output_maximum'][t_idx] #(24)

print("Model setup complete.")
print(f"Number of variables: {len(cqm.variables)}")
print(f"Number of constraints: {len(cqm.constraints)}")

# Solve the CQM
print("Solving...")
sampler = LeapHybridCQMSampler()
solution = sampler.sample_cqm(cqm, time_limit=150)
print("Finished sampling...")

feasible_sampleset = solution.filter(lambda d: d.is_feasible)
if len(feasible_sampleset) > 0:
    best_solution = feasible_sampleset.first
    print("Objective function value:", best_solution.energy)
else:
    print("No feasible solution found.")