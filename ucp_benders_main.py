# Import necessary libraries
import json
from pyomo.environ import *
from master_problem_pyomo import initial_master_problem
from master_problem_pyomo import solving_master_problem
from benders_cuts import generate_benders_cut
from sub_problem import subproblem
import pyomo.environ as pyo
from collections import defaultdict
import random

print("Finished importing")
# Load data
data_file = "rts_gmlc/2020-01-27.json"
print('loading data')
data = json.load(open(data_file, 'r'))

HOURS = 24 #The number of time periods you want
data['time_periods'] = HOURS

# Extract data for generators and time periods
thermal_gens = data['thermal_generators']
renewable_gens = data['renewable_generators']
time_periods = {t+1 : t for t in range(HOURS)}
gen_startup_categories = {g : list(range(0, len(gen['startup']))) for (g, gen) in thermal_gens.items()}
gen_pwl_points = {g : list(range(0, len(gen['piecewise_production']))) for (g, gen) in thermal_gens.items()}



# Iterative Benders Decomposition
max_iterations = 6000
tolerance = 100
iteration = 1
convergence = False
upper_bound = float('inf')
lower_bound = -float('inf')

mumber_of_on_units = 0

'''
counter = 0
for g, gen in thermal_gens.items():
    counter += 1
    if gen['unit_on_t0'] == 0:
        if random.random() < 0.4:
            thermal_gens[g]['unit_on_t0'] = 1 
            thermal_gens[g]['power_output_t0'] = thermal_gens[g]['power_output_minimum']

print(counter)
'''
for g, gen in thermal_gens.items():
    if gen['unit_on_t0'] == 1:
        mumber_of_on_units += 1
print(mumber_of_on_units)
master = initial_master_problem(data, thermal_gens, renewable_gens, time_periods, gen_startup_categories, iteration)



    #mumber_of_on_units = value(gen['unit_on_t0']) + mumber_of_on_units



#Initialyzing loop
print('Solving master problem')
while not convergence and iteration < max_iterations:
    
    #New iteration
    print(f'\nIteration {iteration}')

    #solving master problem
    master_obj_value, master_var_values, beta = solving_master_problem(master, thermal_gens, renewable_gens, time_periods, gen_startup_categories, iteration)

    '''
    constraint_counts = defaultdict(int)
    for constr in master.component_objects(Constraint, active=True):
        for index in constr:
            constraint_counts[constr.name] += 1

    # Print the count for each constraint
    for name, count in constraint_counts.items():
        print(f"Constraint: {name}, Number of constraints: {count}")
    '''

    iteration += 1

    if iteration==2:
        master = initial_master_problem(data, thermal_gens, renewable_gens, time_periods, gen_startup_categories, iteration)
    #setting LB
    lower_bound = master_obj_value
    print(f'Lower bound: {lower_bound}')

    #Solving the sub problem
    dual_values, theta_j = subproblem(data, master_var_values, thermal_gens, renewable_gens, time_periods, gen_pwl_points)

    '''
    for t in time_periods:
        number_of_on_units_after_master = 0
        for g, gen in thermal_gens.items():
            number_of_on_units_after_master += master_var_values["ug"][g, t]
            print(f'number_of_on_units_after_master{t}', number_of_on_units_after_master)
    '''
    
    #Setting UB
    master_cost = master_obj_value - beta
    upper_bound = min(upper_bound, master_cost + theta_j)
    print(f'Upper bound: {upper_bound}')

    print(beta)

    #Finding gap, and setting loop conditions
    gap = upper_bound - lower_bound
    print(f'Optimality gap: {gap}')
    if gap <= tolerance:
        convergence = True
        print('Convergence achieved!')
        break
    if iteration > max_iterations:
        print('max iterations now, stopping algorithm')
        break

    #Generates new benders cuts
    generate_benders_cut(master, dual_values, theta_j, master_var_values)


print("\nFINAL SOLUTION:")
print("The objective value is: ", upper_bound)
print("The number of iterations where: ", iteration, "\n")