# Import necessary libraries
import json
from pyomo.environ import *
from master_problem_pyomo import initial_master_problem
from master_problem_pyomo import solving_master_problem
from benders_cuts import generate_benders_cut
from sub_problem import subproblem
import argparse
import matplotlib.pyplot as plt
import pandas as pd


# Create the parser
parser = argparse.ArgumentParser(description="Run a script with command-line parameters.")

# Add arguments
parser.add_argument("-PENALTY", type=int, help="Set the penalty value", required=True)
parser.add_argument("-BETA_L", type=int, help="Set the beta value", required=True)

# Parse arguments
args = parser.parse_args()

# Access the arguments
PENALTY = args.PENALTY
BETA_L = args.BETA_L

BETA_L =-BETA_L


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
max_iterations = 10
tolerance = 10
iteration = 0
convergence = False
upper_bound = float('inf')
lower_bound = -float('inf')



plotting_values={
    "upper_bound":[],
    "lower_bound":[],
}



master = initial_master_problem(data, thermal_gens, renewable_gens, time_periods, gen_startup_categories, iteration, BETA_L=BETA_L)

#Starting algorithm
print('Solving master problem')
while not convergence and iteration < max_iterations:
    
    #New iteration
    iteration += 1
    print(f'\nIteration {iteration}')

    #solving master problem
    master_obj_value, master_var_values, beta = solving_master_problem(master, thermal_gens, time_periods, gen_startup_categories)

    #if you want to count constraints
    '''
    constraint_counts = defaultdict(int)
    for constr in master.component_objects(Constraint, active=True):
        for index in constr:
            constraint_counts[constr.name] += 1

    # Print the count for each constraint
    for name, count in constraint_counts.items():
        print(f"Constraint: {name}, Number of constraints: {count}")
    '''

    #setting LB
    lower_bound = master_obj_value
    print(f'Lower bound: {lower_bound}')

    #Solving the sub problem
    dual_values, theta_j = subproblem(data, master_var_values, thermal_gens, renewable_gens, time_periods, gen_pwl_points, PENALTY=PENALTY)
    
    #Setting UB
    master_cost = master_obj_value - beta
    upper_bound = min(upper_bound, master_cost + theta_j)
    print(f'Upper bound: {upper_bound}')

    print(f"beta: {beta}")

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

    plotting_values["lower_bound"].append(lower_bound)
    plotting_values["upper_bound"].append(upper_bound)


print("\nFINAL SOLUTION:")
print("The objective value is: ", upper_bound)
print("The number of iterations where: ", iteration, "\n")

plt.plot(plotting_values["lower_bound"], label="lower bound")
plt.plot(plotting_values["upper_bound"], label="upper bound")
plt.yscale("log")
plt.legend()
plt.title(f"Convergence plot. PENALTY:{PENALTY}, BETA_L:{BETA_L}")
plt.savefig(f"plots/plots_penalty_{PENALTY}_beta_l_{BETA_L}.png")
plt.show()


# Save values to CSV
csv_filename = f"csv_files/penalty_{PENALTY}_beta_l_{BETA_L}.csv"
df = pd.DataFrame({
    "lower_bound": plotting_values["lower_bound"],
    "upper_bound": plotting_values["upper_bound"]
})
df.to_csv(csv_filename, index=False)

print(f"CSV file saved to {csv_filename}")