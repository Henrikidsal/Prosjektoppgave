# Import necessary libraries 
from dwave.system import LeapHybridCQMSampler
from dimod import ConstrainedQuadraticModel, Binary, Real
import json
import itertools

print("Finished loading")

# Load data
data_file = "rts_gmlc/2020-01-27.json"
print('Loading data...')
data = json.load(open(data_file, 'r'))

Hours = 12
data["time_periods"] = Hours

# Extract data for generators and time periods
thermal_gens = data['thermal_generators']
renewable_gens = data['renewable_generators']
time_periods = {t + 1: t for t in range(data['time_periods'])}
time_periods = dict(itertools.islice(time_periods.items(), Hours))

gen_startup_categories = {g: list(range(len(gen['startup']))) for g, gen in thermal_gens.items()}
num_pwl_points = 4  # Define the number of piecewise linear points
gen_pwl_points = {
    g: list(range(min(num_pwl_points, len(gen['piecewise_production']))))
    for g, gen in thermal_gens.items()
}

print('Building model...')
# Initialize the Constrained Quadratic Model (CQM)
cqm = ConstrainedQuadraticModel()

# Define Variables

# Continuous Variables for thermal generators
cg = {}  # Generation cost
pg = {}  # Power generation above minimum
rg = {}  # Reserves
# Binary Variables for thermal generators
ug = {}  # Unit on/off
vg = {}  # Unit startup
wg = {}  # Unit shutdown
# Binary Variables for startup categories
dg = {}
# Continuous Variables for piecewise linear generation
lg = {}
# Variables for renewable generators
pw = {}  # Power generation

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

# Demand Constraint (2)
for t, t_idx in time_periods.items():
    demand_expr = sum(
        pg[g, t] + thermal_gens[g]['power_output_minimum'] * ug[g, t]
        for g in thermal_gens
    )
    demand_expr += sum(pw[w, t] for w in renewable_gens)
    cqm.add_constraint(
        demand_expr == data['demand'][t_idx], label=f'demand_constraint_{t}'
    )

# Reserves Constraint (3)
for t, t_idx in time_periods.items():
    reserve_expr = sum(rg[g, t] for g in thermal_gens)
    cqm.add_constraint(
        reserve_expr >= data['reserves'][t_idx], label=f'reserve_constraint_{t}'
    )

# Initial Conditions and Ramp Limits Constraints
for g, gen in thermal_gens.items():
    # Uptime and Downtime Constraints at t=0
    if gen['unit_on_t0'] == 1:
        if gen['time_up_minimum'] - gen['time_up_t0'] >= 1:
            UT = min(
                gen['time_up_minimum'] - gen['time_up_t0'], data['time_periods']
            )
            if UT >= 1:
                uptime_expr = sum(
                    (ug[g, t] - 1) for t in range(1, UT + 1)
                )
                cqm.add_constraint(
                    uptime_expr == 0, label=f'uptime_{g}_t0'
                )
    elif gen['unit_on_t0'] == 0:
        if gen['time_down_minimum'] - gen['time_down_t0'] >= 1:
            DT = min(
                gen['time_down_minimum'] - gen['time_down_t0'],
                data['time_periods'],
            )
            if DT >= 1:
                downtime_expr = sum(
                    ug[g, t] for t in range(1, DT + 1)
                )
                cqm.add_constraint(
                    downtime_expr == 0, label=f'downtime_{g}_t0'
                )
    else:
        raise Exception(
            f'Invalid unit_on_t0 for generator {g}, unit_on_t0={gen["unit_on_t0"]}'
        )

    # Logical Constraints at t=1 (6)
    expr = ug[g, 1] - gen['unit_on_t0'] - vg[g, 1] + wg[g, 1]
    cqm.add_constraint(expr == 0, label=f'logical_{g}_t0')

    # Startup Allowed Constraints at t=0 (7)
    startup_terms = []
    for s in gen_startup_categories[g][:-1]:  # All but last
        lag_s1 = gen['startup'][s + 1]['lag']
        lag_s = gen['startup'][s]['lag']
        start = max(1, lag_s1 - gen['time_down_t0'] + 1)
        end = min(lag_s1 - 1, data['time_periods'])
        if start <= end:
            for t in range(start, end + 1):
                startup_terms.append(dg[g, s, t])
    if startup_terms:
        startup_expr = sum(startup_terms)
        cqm.add_constraint(
            startup_expr == 0, label=f'startup_allowed_{g}_t0'
        )

    # Ramp Up Limit at t=1 (8)
    ramp_up_expr = (
        pg[g, 1] + rg[g, 1] - gen['unit_on_t0'] *
        (gen['power_output_t0'] - gen['power_output_minimum'])
    )
    cqm.add_constraint(
        ramp_up_expr <= gen['ramp_up_limit'], label=f'ramp_up_{g}_t0'
    )

    # Ramp Down Limit at t=1 (9)
    ramp_down_expr = (
        gen['unit_on_t0'] * (gen['power_output_t0'] - gen['power_output_minimum'])
        - pg[g, 1]
    )
    cqm.add_constraint(
        ramp_down_expr <= gen['ramp_down_limit'], label=f'ramp_down_{g}_t0'
    )

    # Shutdown Constraint at t=1 (10)
    Pmax = gen['power_output_maximum'] - gen['power_output_minimum']
    ramp_sd_limit = max(
        (gen['power_output_maximum'] - gen['ramp_shutdown_limit']), 0
    )
    shutdown_expr = (
        gen['unit_on_t0'] * (gen['power_output_t0'] - gen['power_output_minimum'])
        - (gen['unit_on_t0'] * Pmax - ramp_sd_limit * wg[g, 1])
    )
    cqm.add_constraint(
        shutdown_expr <= 0, label=f'shutdown_{g}_t0'
    )

# Must-run Constraints (11) and other operational constraints
for g, gen in thermal_gens.items():
    UT_g = min(gen['time_up_minimum'], data['time_periods'])
    DT_g = min(gen['time_down_minimum'], data['time_periods'])
    for t in time_periods:
        # Must-run Constraint (11)
        cqm.add_constraint(
            ug[g, t] >= gen['must_run'], label=f'must_run_{g}_{t}'
        )

        if t > 1:
            # Logical Constraint (12)
            expr = ug[g, t] - ug[g, t - 1] - vg[g, t] + wg[g, t]
            cqm.add_constraint(expr == 0, label=f'logical_{g}_{t}')

            # Ramp Up Limit (19)
            ramp_up_expr = pg[g, t] + rg[g, t] - pg[g, t - 1]
            cqm.add_constraint(
                ramp_up_expr <= gen['ramp_up_limit'], label=f'ramp_up_{g}_{t}'
            )

            # Ramp Down Limit (20)
            ramp_down_expr = pg[g, t - 1] - pg[g, t]
            cqm.add_constraint(
                ramp_down_expr <= gen['ramp_down_limit'], label=f'ramp_down_{g}_{t}'
            )

    # Uptime Constraints (13)
    if UT_g > 0:
        for t in range(min(UT_g, data['time_periods']), data['time_periods'] + 1):
            uptime_expr = sum(vg[g, tp] for tp in range(t - UT_g + 1, t + 1)) - ug[g, t]
            cqm.add_constraint(uptime_expr <= 0, label=f'uptime_{g}_{t}')

    # Downtime Constraints (14)
    if DT_g > 0:
        for t in range(min(DT_g, data['time_periods']), data['time_periods'] + 1):
            downtime_expr = sum(wg[g, tp] for tp in range(t - DT_g + 1, t + 1)) - (1 - ug[g, t])
            cqm.add_constraint(downtime_expr <= 0, label=f'downtime_{g}_{t}')

    # Startup Allowed Constraints (15)
    for s in gen_startup_categories[g][:-1]:  # All but last
        lag_s = gen['startup'][s]['lag']
        lag_s1 = gen['startup'][s + 1]['lag']
        for t in range(lag_s1, data['time_periods'] + 1):
            startup_allowed_terms = []
            for i in range(lag_s, lag_s1):
                time_index = t - i
                if time_index >= 1 and time_index in time_periods:
                    startup_allowed_terms.append(wg[g, time_index])
            if startup_allowed_terms:
                startup_allowed_expr = (
                    dg[g, s, t] - sum(startup_allowed_terms)
                )
                cqm.add_constraint(
                    startup_allowed_expr <= 0, label=f'startup_allowed_{g}_{s}_{t}'
                )

    for t in time_periods:
        # Startup Selection Constraint (16)
        startup_select_expr = vg[g, t] - sum(
            dg[g, s, t] for s in gen_startup_categories[g]
        )
        cqm.add_constraint(
            startup_select_expr == 0, label=f'startup_select_{g}_{t}'
        )

        # Generation Limits Constraints (17) and (18)
        Pmax = gen['power_output_maximum'] - gen['power_output_minimum']
        ramp_su_limit = max(
            (gen['power_output_maximum'] - gen['ramp_startup_limit']), 0
        )
        ramp_sd_limit = max(
            (gen['power_output_maximum'] - gen['ramp_shutdown_limit']), 0
        )

        # Generation Limit Constraint 1 (17)
        gen_limit_expr1 = (
            pg[g, t] + rg[g, t] - (Pmax * ug[g, t] - ramp_su_limit * vg[g, t])
        )
        cqm.add_constraint(
            gen_limit_expr1 <= 0, label=f'gen_limit1_{g}_{t}'
        )

        # Generation Limit Constraint 2 (18)
        if t < data['time_periods']:
            next_t = t + 1
            gen_limit_expr2 = (
                pg[g, t] + rg[g, t] - (Pmax * ug[g, t] - ramp_sd_limit * wg[g, next_t])
            )
            cqm.add_constraint(
                gen_limit_expr2 <= 0, label=f'gen_limit2_{g}_{t}'
            )

        # Power Selection Constraint (21)
        piece_mw0 = gen['piecewise_production'][0]['mw']
        power_select_expr = pg[g, t] - sum(
            (gen['piecewise_production'][l]['mw'] - piece_mw0) * lg[g, l, t]
            for l in gen_pwl_points[g]
        )
        cqm.add_constraint(
            power_select_expr == 0, label=f'power_select_{g}_{t}'
        )

        # Cost Selection Constraint (22)
        piece_cost0 = gen['piecewise_production'][0]['cost']
        cost_select_expr = cg[g, t] - sum(
            (gen['piecewise_production'][l]['cost'] - piece_cost0) * lg[g, l, t]
            for l in gen_pwl_points[g]
        )
        cqm.add_constraint(
            cost_select_expr == 0, label=f'cost_select_{g}_{t}'
        )

        # On Selection Constraint (23)
        on_select_expr = ug[g, t] - sum(
            lg[g, l, t] for l in gen_pwl_points[g]
        )
        cqm.add_constraint(
            on_select_expr == 0, label=f'on_select_{g}_{t}'
        )

# Renewable Generators Constraints (24)
for w, gen in renewable_gens.items():
    for t, t_idx in time_periods.items():
        # Power output limits
        pw_min = gen['power_output_minimum'][t_idx]
        pw_max = gen['power_output_maximum'][t_idx]
        cqm.add_constraint(
            pw[w, t] >= pw_min, label=f'renewable_min_{w}_{t}'
        )
        cqm.add_constraint(
            pw[w, t] <= pw_max, label=f'renewable_max_{w}_{t}'
        )

print("Model setup complete.")
print(f"Number of variables: {len(cqm.variables)}")
print(f"Number of constraints: {len(cqm.constraints)}")
# Solve the CQM
print("Solving...")
sampler = LeapHybridCQMSampler()
solution = sampler.sample_cqm(cqm, time_limit=400)
print("Finished sampling...")

# Process and display solution
feasible_sampleset = solution.filter(lambda d: d.is_feasible)
if len(feasible_sampleset) > 0:
    best_solution = feasible_sampleset.first
    print("Feasible solution found with energy:", best_solution.energy)
    # Extract the variable values
    # sample = best_solution.sample
    # Print the variable values
    # for var, val in sample.items():
    #     print(f"{var}: {val}")
else:
    print("No feasible solution found.")
