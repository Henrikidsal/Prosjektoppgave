# Import necessary libraries
import json
from pyomo.environ import *
from master_problem_pyomo import Master_problem_pyomo
from sub_problem import subproblem
from benders_cuts import generate_benders_cut

print("Finished importing")
# Load data
data_file = "rts_gmlc/2020-01-27.json"
print('loading data')
data = json.load(open(data_file, 'r'))

HOURS = 48 #The number of time periods you want
data['time_periods'] = HOURS

# Extract data for generators and time periods
thermal_gens = data['thermal_generators']
renewable_gens = data['renewable_generators']
time_periods = {t+1 : t for t in range(HOURS)}
gen_startup_categories = {g : list(range(0, len(gen['startup']))) for (g, gen) in thermal_gens.items()}
gen_pwl_points = {g : list(range(0, len(gen['piecewise_production']))) for (g, gen) in thermal_gens.items()}

# Algorithm for solving:
max_iterations = 6000
tolerance = 100
iteration = 0
convergence = False
upper_bound = float('inf')
lower_bound = -float('inf')



master,  =initialize_master()

while not convergence and iteration < max_iterations:
        
    iteration += 1
    print(f'\nIteration {iteration}')

    #solve master problem
    master_obj_value, master_var_values, beta, master = Master_problem_pyomo(data, thermal_gens, time_periods, gen_startup_categories)
    
    #set LB
    lower_bound = master_obj_value
    print(f'Lower bound: {lower_bound}')

    #Solves sub problem
    dual_values, sub_cost = subproblem(data, master_var_values, thermal_gens, renewable_gens, time_periods, gen_pwl_points)

    #Setting UB
    master_cost = master_obj_value - beta
    upper_bound = min(upper_bound, master_cost + sub_cost)
    print(f'Upper bound: {upper_bound}')

    #Finding gap, and setting loop conditions
    gap = upper_bound - lower_bound
    print(f'Optimality gap: {gap}')
    if gap <= tolerance:
        convergence = True
        print('Convergence achieved!')
        break
    if iteration >= max_iterations:
        print('Maximum iterations reached without convergence.')
        break

    #Generate new benders cuts
    generate_benders_cut(master, dual_values, sub_cost, master_var_values)


print("\nFINAL SOLUTION:")
print("The objective value is: ", upper_bound)
print(f'Optimality gap: {gap}')
print("The number of iterations where: ", iteration, "\n")