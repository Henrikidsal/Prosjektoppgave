from pyomo.environ import *
import json
from pyomo.opt import SolverFactory
import copy
import random
import matplotlib.pyplot as plt
import numpy as np
# Load data
data_file = "rts_gmlc/2020-11-25.json"
print('loading data')
data = json.load(open(data_file, 'r'))

#INPUTS:
random.seed(19)

# How many HOURS do you want in the problem?
HOURS = 48

# How much do you want to reduce generator capasity and demand?
reduction_percentage = 0.8

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

    # Scale reserves proportionally to the scaled demand
    total_initial_demand = sum(demand)
    total_scaled_demand = sum(scaled_demand)
    reserve_scale_factor = total_scaled_demand / total_initial_demand if total_initial_demand > 0 else 1

    scaled_reserves = [r * reserve_scale_factor for r in reserves]

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
m = ConcreteModel()

m.cg = Var(thermal_gens.keys(), time_periods.keys())
m.pg = Var(thermal_gens.keys(), time_periods.keys(), within=NonNegativeReals)  
m.rg = Var(thermal_gens.keys(), time_periods.keys(), within=NonNegativeReals)  
m.pw = Var(renewable_gens.keys(), time_periods.keys(), within=NonNegativeReals)
m.ug = Var(thermal_gens.keys(), time_periods.keys(), within=Binary) 
m.vg = Var(thermal_gens.keys(), time_periods.keys(), within=Binary) 
m.wg = Var(thermal_gens.keys(), time_periods.keys(), within=Binary) 
m.dg = Var(((g,s,t) for g in thermal_gens for s in gen_startup_categories[g] for t in time_periods), within=Binary) ##
m.lg = Var(((g,l,t) for g in thermal_gens for l in gen_pwl_points[g] for t in time_periods), within=UnitInterval) ##

m.obj = Objective(expr=sum(
                sum(
                    m.cg[g,t] + gen['piecewise_production'][0]['cost']*m.ug[g,t]
                    + sum( gen_startup['cost']*m.dg[g,s,t] for (s, gen_startup) in enumerate(gen['startup']))
                for t in time_periods)
            for g, gen in thermal_gens.items() )
            )#(1)

m.demand = Constraint(time_periods.keys())
m.reserves = Constraint(time_periods.keys())
for t,t_idx in time_periods.items():
    m.demand[t] = sum( m.pg[g,t]+gen['power_output_minimum']*m.ug[g,t] for (g, gen) in thermal_gens.items() ) + sum( m.pw[w,t] for w in renewable_gens ) == demand[t_idx] #(2)
    m.reserves[t] = sum( m.rg[g,t] for g in thermal_gens ) >= reserves[t_idx] #(3)

m.uptimet0 = Constraint(thermal_gens.keys())
m.downtimet0 = Constraint(thermal_gens.keys())
m.logicalt0 = Constraint(thermal_gens.keys())
m.startupt0 = Constraint(thermal_gens.keys())

m.rampupt0 = Constraint(thermal_gens.keys())
m.rampdownt0 = Constraint(thermal_gens.keys())
m.shutdownt0 = Constraint(thermal_gens.keys())

for g, gen in thermal_gens.items():
    if gen['unit_on_t0'] == 1:
        if gen['time_up_minimum'] - gen['time_up_t0'] >= 1:
            m.uptimet0[g] = sum( (m.ug[g,t] - 1) for t in range(1, min(gen['time_up_minimum'] - gen['time_up_t0'], HOURS)+1)) == 0 #(4)
    elif gen['unit_on_t0'] == 0:
        if gen['time_down_minimum'] - gen['time_down_t0'] >= 1:
            m.downtimet0[g] = sum( m.ug[g,t] for t in range(1, min(gen['time_down_minimum'] - gen['time_down_t0'], HOURS)+1)) == 0 #(5)
    else:
        raise Exception('Invalid unit_on_t0 for generator {}, unit_on_t0={}'.format(g, gen['unit_on_t0']))

    m.logicalt0[g] = m.ug[g,1] - gen['unit_on_t0'] == m.vg[g,1] - m.wg[g,1] #(6)

    startup_expr = sum( 
                        sum( m.dg[g,s,t] 
                                for t in range(
                                                max(1, gen['startup'][s+1]['lag']-gen['time_down_t0']+1),
                                                min(gen['startup'][s+1]['lag']-1, HOURS)+1
                                              )
                            )
                       for s,_ in enumerate(gen['startup'][:-1])) ## all but last
    if isinstance(startup_expr, int):
        pass
    else:
        m.startupt0[g] = startup_expr == 0 #(7)

    m.rampupt0[g] = m.pg[g,1] + m.rg[g,1] - gen['unit_on_t0']*(gen['power_output_t0']-gen['power_output_minimum']) <= gen['ramp_up_limit'] #(8)

    m.rampdownt0[g] = gen['unit_on_t0']*(gen['power_output_t0']-gen['power_output_minimum']) - m.pg[g,1] <= gen['ramp_down_limit'] #(9)

    shutdown_constr = gen['unit_on_t0']*(gen['power_output_t0']-gen['power_output_minimum']) <= gen['unit_on_t0']*(gen['power_output_maximum'] - gen['power_output_minimum']) - max((gen['power_output_maximum'] - gen['ramp_shutdown_limit']),0)*m.wg[g,1] #(10)

    if isinstance(shutdown_constr, bool):
        pass
    else:
        m.shutdownt0[g] = shutdown_constr

m.mustrun = Constraint(thermal_gens.keys(), time_periods.keys())
m.logical = Constraint(thermal_gens.keys(), time_periods.keys())
m.uptime = Constraint(thermal_gens.keys(), time_periods.keys())
m.downtime = Constraint(thermal_gens.keys(), time_periods.keys())
m.startup_select = Constraint(thermal_gens.keys(), time_periods.keys())
m.gen_limit1 = Constraint(thermal_gens.keys(), time_periods.keys())
m.gen_limit2 = Constraint(thermal_gens.keys(), time_periods.keys())
m.ramp_up = Constraint(thermal_gens.keys(), time_periods.keys())
m.ramp_down = Constraint(thermal_gens.keys(), time_periods.keys())
m.power_select = Constraint(thermal_gens.keys(), time_periods.keys())
m.cost_select = Constraint(thermal_gens.keys(), time_periods.keys())
m.on_select = Constraint(thermal_gens.keys(), time_periods.keys())

for g, gen in thermal_gens.items():
    for t in time_periods:
        m.mustrun[g,t] = m.ug[g,t] >= gen['must_run'] #(11)

        if t > 1:
            m.logical[g,t] = m.ug[g,t] - m.ug[g,t-1] == m.vg[g,t] - m.wg[g,t] #(12)

        UT = min(gen['time_up_minimum'], HOURS)
        if t >= UT:
            m.uptime[g,t] = sum(m.vg[g,t] for t in range(t-UT+1, t+1)) <= m.ug[g,t] #(13)
        DT = min(gen['time_down_minimum'], HOURS)
        if t >= DT:
            m.downtime[g,t] = sum(m.wg[g,t] for t in range(t-DT+1, t+1)) <= 1-m.ug[g,t] #(14)
        m.startup_select[g,t] = m.vg[g,t] == sum(m.dg[g,s,t] for s,_ in enumerate(gen['startup'])) #(16)

        m.gen_limit1[g,t] = m.pg[g,t]+m.rg[g,t] <= (gen['power_output_maximum'] - gen['power_output_minimum'])*m.ug[g,t] - max((gen['power_output_maximum'] - gen['ramp_startup_limit']),0)*m.vg[g,t] #(17)

        if t < len(time_periods): 
            m.gen_limit2[g,t] = m.pg[g,t]+m.rg[g,t] <= (gen['power_output_maximum'] - gen['power_output_minimum'])*m.ug[g,t] - max((gen['power_output_maximum'] - gen['ramp_shutdown_limit']),0)*m.wg[g,t+1] #(18)

        if t > 1:
            m.ramp_up[g,t] = m.pg[g,t]+m.rg[g,t] - m.pg[g,t-1] <= gen['ramp_up_limit'] #(19)
            m.ramp_down[g,t] = m.pg[g,t-1] - m.pg[g,t] <= gen['ramp_down_limit'] #(20

        piece_mw1 = gen['piecewise_production'][0]['mw']
        piece_cost1 = gen['piecewise_production'][0]['cost']
        m.power_select[g,t] = m.pg[g,t] == sum( (piece['mw'] - piece_mw1)*m.lg[g,l,t] for l,piece in enumerate(gen['piecewise_production'])) #(21)
        m.cost_select[g,t] = m.cg[g,t] == sum( (piece['cost'] - piece_cost1)*m.lg[g,l,t] for l,piece in enumerate(gen['piecewise_production'])) #(22)
        m.on_select[g,t] = m.ug[g,t] == sum(m.lg[g,l,t] for l,_ in enumerate(gen['piecewise_production'])) #(23)


m.startup_allowed = Constraint(m.dg.index_set())
for g, gen in thermal_gens.items():
    for s,_ in enumerate(gen['startup'][:-1]): ## all but last
        for t in time_periods:
            if t >= gen['startup'][s+1]['lag']:
                m.startup_allowed[g,s,t] = m.dg[g,s,t] <= sum(m.wg[g,t-i] for i in range(gen['startup'][s]['lag'], gen['startup'][s+1]['lag'])) #(15)

for w, gen in renewable_gens.items():
    for t, t_idx in time_periods.items():
        m.pw[w,t].setlb(gen['power_output_minimum'][t_idx]) #(24)
        m.pw[w,t].setub(gen['power_output_maximum'][t_idx]) #(24)

print("model setup complete")

gurobi = SolverFactory('gurobi')

print("solving")
result = gurobi.solve(m, options={'MIPGap':0.0}, tee=True)

num_variables = sum(1 for _ in m.component_data_objects(Var, active=True))
print(f"Number of variables: {num_variables}")

num_constraints = sum(len(constraint) for constraint in m.component_objects(Constraint, active=True))
print(f"Number of constraints: {num_constraints}")

m.solutions.load_from(result)
print(f"Objective function value: {value(m.obj)}")


