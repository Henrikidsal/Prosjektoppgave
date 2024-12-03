# Import necessary libraries
import json
from pyomo.environ import *
from master_problem_pyomo import Master_problem_pyomo
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




lower_bound, upper_bound, iterations = Master_problem_pyomo(data, thermal_gens, renewable_gens, time_periods, gen_startup_categories, gen_pwl_points)
print("\nFINAL SOLUTION:")
print("The objective value is: ", upper_bound)
print("The number of iterations where: ", iterations, "\n")