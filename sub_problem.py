# Import necessary libraries
import numpy as np
from pyomo.environ import *
from pyomo.opt import SolverFactory
import pyomo.environ as pyo
from pyomo.environ import Suffix
from pyomo.opt import TerminationCondition

# Subproblem
def subproblem(data, master_solution, thermal_gens, renewable_gens, time_periods, gen_pwl_points):

    m = ConcreteModel()
    print("creating sub problem")

    #sub problem variables
    m.cg = Var(thermal_gens.keys(), time_periods.keys())
    m.pg = Var(thermal_gens.keys(), time_periods.keys(), within=NonNegativeReals)  
    m.rg = Var(thermal_gens.keys(), time_periods.keys(), within=NonNegativeReals)  
    m.pw = Var(renewable_gens.keys(), time_periods.keys(), within=NonNegativeReals)
    m.lg = Var(((g,l,t) for g in thermal_gens for l in gen_pwl_points[g] for t in time_periods), within=UnitInterval) ##
    
    # Parameters (Fixed from master problem)
    m.ug = Param(thermal_gens.keys(), time_periods.keys(), initialize=master_solution['ug'], within=Binary)
    m.vg = Param(thermal_gens.keys(), time_periods.keys(), initialize=master_solution['vg'], within=Binary)
    m.wg = Param(thermal_gens.keys(), time_periods.keys(), initialize=master_solution['wg'], within=Binary)

    # Objective Function
    m.obj = Objective(expr=sum(sum(m.cg[g, t] for t in time_periods) for g in thermal_gens))

    # Constraints
    def demand_rule(m, t):
        t_idx = time_periods[t]
        return sum( m.pg[g,t]+gen['power_output_minimum']*m.ug[g,t] for (g, gen) in thermal_gens.items() ) + sum( m.pw[w,t] for w in renewable_gens ) == data['demand'][t_idx] #(2)
    m.demand = Constraint(time_periods.keys(), rule=demand_rule, doc='demand')

    def reserves_rule(m, t):
        t_idx = time_periods[t]
        return sum( m.rg[g,t] for g in thermal_gens ) >= data['reserves'][t_idx] #(3)
    m.reserves = Constraint(time_periods.keys(), rule=reserves_rule, doc='reserves')

    def rampupt0_rule(m, g):
        gen = thermal_gens[g]
        return (m.pg[g,1] + m.rg[g,1] - gen['unit_on_t0']*(gen['power_output_t0']-gen['power_output_minimum']) <= gen['ramp_up_limit']) #(8)
    m.rampupt0 = Constraint(thermal_gens.keys(), rule=rampupt0_rule, doc='rampupt0')

    def rampdownt0_rule(m, g):
        gen = thermal_gens[g]
        return (gen['unit_on_t0']*(gen['power_output_t0']-gen['power_output_minimum']) - m.pg[g,1] <= gen['ramp_down_limit']) #(9)
    m.rampdownt0 = Constraint(thermal_gens.keys(), rule=rampdownt0_rule, doc='rampdownt0')

    def gen_limit1_rule(m, g, t):
        gen = thermal_gens[g]
        return (m.pg[g,t]+m.rg[g,t] <= (gen['power_output_maximum'] - gen['power_output_minimum'])*m.ug[g,t] - max((gen['power_output_maximum'] - gen['ramp_startup_limit']),0)*m.vg[g,t]) #(17)
    m.gen_limit1 = Constraint(thermal_gens.keys(), time_periods.keys(), rule=gen_limit1_rule, doc='gen_limit1')

    def gen_limit2_rule(m, g, t):
        gen = thermal_gens[g]
        if t < len(time_periods):
            return (m.pg[g,t]+m.rg[g,t] <= (gen['power_output_maximum'] - gen['power_output_minimum'])*m.ug[g,t] - max((gen['power_output_maximum'] - gen['ramp_shutdown_limit']),0)*m.wg[g,t+1])
        else:
            return Constraint.Skip
    m.gen_limit2 = Constraint(thermal_gens.keys(), time_periods.keys(), rule=gen_limit2_rule, doc='gen_limit2')

    def ramp_up_rule(m, g, t):
        gen = thermal_gens[g]
        if t > 1:
            return m.pg[g,t]+m.rg[g,t] - m.pg[g,t-1] <= gen['ramp_up_limit'] #(19)
        else:
            return Constraint.Skip
    m.ramp_up = Constraint(thermal_gens.keys(), time_periods.keys(), rule=ramp_up_rule, doc='ramp-up')

    def ramp_down_rule(m, g, t):
        gen = thermal_gens[g]
        if t > 1:
            return m.pg[g,t-1] - m.pg[g,t] <= gen['ramp_down_limit'] #(20)
        else:
            return Constraint.Skip
    m.ramp_down = Constraint(thermal_gens.keys(), time_periods.keys(), rule=ramp_down_rule, doc='ramp-down')

    # Piecewise Linear Constraints
    def power_select_rule(m, g, t):
        gen = thermal_gens[g]
        piece_mw1 = gen['piecewise_production'][0]['mw']
        return m.pg[g,t] == sum( (piece['mw'] - piece_mw1)*m.lg[g,l,t] for l,piece in enumerate(gen['piecewise_production'])) #(21)
    m.power_select = Constraint(thermal_gens.keys(), time_periods.keys(), rule=power_select_rule, doc='power_select')

    def cost_select_rule(m, g, t):
        gen = thermal_gens[g]
        piece_cost1 = gen['piecewise_production'][0]['cost']
        return m.cg[g,t] == sum( (piece['cost'] - piece_cost1)*m.lg[g,l,t] for l,piece in enumerate(gen['piecewise_production'])) #(22)
    m.cost_select = Constraint(thermal_gens.keys(), time_periods.keys(), rule=cost_select_rule, doc='cost_select')

    def on_select_rule(m, g, t):
        return sum(m.lg[g,l,t] for l in gen_pwl_points[g]) == m.ug[g,t]
    m.on_select = Constraint(thermal_gens.keys(), time_periods.keys(), rule=on_select_rule, doc='on_select')


    for w, gen in renewable_gens.items():
        for t, t_idx in time_periods.items():
            m.pw[w,t].setlb(gen['power_output_minimum'][t_idx]) #(24)
            m.pw[w,t].setub(gen['power_output_maximum'][t_idx]) #(24)

    #Collecting dual variables
    m.dual = Suffix(direction=Suffix.IMPORT)

    print('Solving subproblem')
    solver = SolverFactory('gurobi')
    solver.options['OutputFlag'] = 0
    solver.options['DualReductions'] = 0
    results = solver.solve(m, tee=False)

    solver_status = True
    dual_values = {}
    if results.solver.termination_condition == TerminationCondition.infeasible:
        print("Subproblem is infeasible.")
        solver_status = False
    elif results.solver.termination_condition == TerminationCondition.optimal:
        print("Subproblem is feasible.")
        solver_status = True
    else:
        print("Solver did not return optimality or infeasibility.")
        solver_status = None

    if solver_status:
        try:
            for c in [m.demand, m.reserves, m.rampupt0, m.rampdownt0, m.gen_limit1, m.gen_limit2, m.ramp_up, m.ramp_down, m.power_select, m.cost_select, m.on_select]:
                for index in c:
                    con = c[index]
                    if con in m.dual:
                        dual_values[index] = m.dual[con]
                    else:
                        dual_values[index] = None
        except KeyError as e:
            print(f"Error accessing duals: {e}")
            dual_values = None
    else:
        print("Solving feasibility dual problem to obtain dual variables for feasibility cut.")
        dual_values = solve_feasibility_dual(data, master_solution, thermal_gens, renewable_gens, time_periods, gen_pwl_points)

    return m, dual_values, solver_status

'''
def solve_feasibility_dual(data, master_solution, thermal_gens, renewable_gens, time_periods, gen_pwl_points):
    # Build the primal subproblem model (similar to your subproblem but without the objective)
    m = ConcreteModel()

    # Variables
    m.cg = Var(thermal_gens.keys(), time_periods.keys())
    m.pg = Var(thermal_gens.keys(), time_periods.keys(), within=NonNegativeReals)  
    m.rg = Var(thermal_gens.keys(), time_periods.keys(), within=NonNegativeReals)  
    m.pw = Var(renewable_gens.keys(), time_periods.keys(), within=NonNegativeReals)
    m.lg = Var(((g,l,t) for g in thermal_gens for l in gen_pwl_points[g] for t in time_periods), within=UnitInterval)

    # Parameters (Fixed from master problem)
    m.ug = Var(thermal_gens.keys(), time_periods.keys(), within=NonNegativeReals)
    m.vg = Var(thermal_gens.keys(), time_periods.keys(), within=NonNegativeReals)
    m.wg = Var(thermal_gens.keys(), time_periods.keys(), within=NonNegativeReals)

    for g in thermal_gens.keys():
        for t in time_periods.keys():
            m.ug[g, t].fix(master_solution['ug'][(g, t)])
            m.vg[g, t].fix(master_solution['vg'][(g, t)])
            m.wg[g, t].fix(master_solution['wg'][(g, t)])

    def demand_rule(m, t):
        t_idx = time_periods[t]
        return sum(m.pg[g, t] + thermal_gens[g]['power_output_minimum'] * m.ug[g, t] for g in thermal_gens) \
               + sum(m.pw[w, t] for w in renewable_gens) == data['demand'][t_idx]
    m.demand = Constraint(time_periods.keys(), rule=demand_rule, doc='demand')

    def reserves_rule(m, t):
        t_idx = time_periods[t]
        return sum(m.rg[g, t] for g in thermal_gens) >= data['reserves'][t_idx]
    m.reserves = Constraint(time_periods.keys(), rule=reserves_rule, doc='reserves')

    def rampupt0_rule(m, g):
        gen = thermal_gens[g]
        return (m.pg[g, 1] + m.rg[g, 1] - gen['unit_on_t0'] * (gen['power_output_t0'] - gen['power_output_minimum']) <= gen['ramp_up_limit'])
    m.rampupt0 = Constraint(thermal_gens.keys(), rule=rampupt0_rule, doc='rampupt0')

    def rampdownt0_rule(m, g):
        gen = thermal_gens[g]
        return (gen['unit_on_t0'] * (gen['power_output_t0'] - gen['power_output_minimum']) - m.pg[g, 1] <= gen['ramp_down_limit'])
    m.rampdownt0 = Constraint(thermal_gens.keys(), rule=rampdownt0_rule, doc='rampdownt0')

    def gen_limit1_rule(m, g, t):
        gen = thermal_gens[g]
        return (m.pg[g, t] + m.rg[g, t] <= (gen['power_output_maximum'] - gen['power_output_minimum']) * m.ug[g, t] -
                max((gen['power_output_maximum'] - gen['ramp_startup_limit']), 0) * m.vg[g, t])
    m.gen_limit1 = Constraint(thermal_gens.keys(), time_periods.keys(), rule=gen_limit1_rule, doc='gen_limit1')

    def gen_limit2_rule(m, g, t):
        gen = thermal_gens[g]
        if t < len(time_periods):
            return (m.pg[g, t] + m.rg[g, t] <= (gen['power_output_maximum'] - gen['power_output_minimum']) * m.ug[g, t] -
                    max((gen['power_output_maximum'] - gen['ramp_shutdown_limit']), 0) * m.wg[g, t + 1])
        else:
            return Constraint.Skip
    m.gen_limit2 = Constraint(thermal_gens.keys(), time_periods.keys(), rule=gen_limit2_rule, doc='gen_limit2')

    def ramp_up_rule(m, g, t):
        gen = thermal_gens[g]
        if t > 1:
            return m.pg[g, t] + m.rg[g, t] - m.pg[g, t - 1] <= gen['ramp_up_limit']
        else:
            return Constraint.Skip
    m.ramp_up = Constraint(thermal_gens.keys(), time_periods.keys(), rule=ramp_up_rule, doc='ramp-up')

    def ramp_down_rule(m, g, t):
        gen = thermal_gens[g]
        if t > 1:
            return m.pg[g, t - 1] - m.pg[g, t] <= gen['ramp_down_limit']
        else:
            return Constraint.Skip
    m.ramp_down = Constraint(thermal_gens.keys(), time_periods.keys(), rule=ramp_down_rule, doc='ramp-down')

    # Piecewise Linear Constraints
    def power_select_rule(m, g, t):
        gen = thermal_gens[g]
        piece_mw1 = gen['piecewise_production'][0]['mw']
        return m.pg[g, t] == sum((piece['mw'] - piece_mw1) * m.lg[g, l, t] for l, piece in enumerate(gen['piecewise_production']))
    m.power_select = Constraint(thermal_gens.keys(), time_periods.keys(), rule=power_select_rule, doc='power_select')

    def cost_select_rule(m, g, t):
        gen = thermal_gens[g]
        piece_cost1 = gen['piecewise_production'][0]['cost']
        return m.cg[g, t] == sum((piece['cost'] - piece_cost1) * m.lg[g, l, t] for l, piece in enumerate(gen['piecewise_production']))
    m.cost_select = Constraint(thermal_gens.keys(), time_periods.keys(), rule=cost_select_rule, doc='cost_select')

    def on_select_rule(m, g, t):
        gen = thermal_gens[g]
        return m.ug[g, t] == sum(m.lg[g, l, t] for l, _ in enumerate(gen['piecewise_production']))
    m.on_select = Constraint(thermal_gens.keys(), time_periods.keys(), rule=on_select_rule, doc='on_select')

    # Bounds for renewable generators
    for w, gen in renewable_gens.items():
        for t, t_idx in time_periods.items():
            m.pw[w, t].setlb(gen['power_output_minimum'][t_idx])
            m.pw[w, t].setub(gen['power_output_maximum'][t_idx])

    m.obj = Objective(expr=0)

    # Dual Suffix
    m.dual = Suffix(direction=Suffix.IMPORT_EXPORT)

    # Apply the duality transformation
    dual = TransformationFactory('duality.linear_dual').create_using(m)

    # Solve the dual model
    print('Solving feasibility dual problem')
    solver = SolverFactory('gurobi')
    solver.options['OutputFlag'] = 0
    results = solver.solve(dual, tee=False)

    feasibility_duals = {}
    if results.solver.termination_condition == TerminationCondition.optimal:
        print("Feasibility dual problem solved optimally.")
        try:
            # Extract dual variables from the dual model
            for var in dual.component_objects(Var, active=True):
                varname = var.getname()
                for index in var:
                    feasibility_duals[(varname, index)] = var[index].value
        except Exception as e:
            print(f"Error accessing dual variables: {e}")
            feasibility_duals = None
    else:
        print("Feasibility dual problem did not solve optimally or is infeasible.")
        feasibility_duals = None

    return feasibility_duals
'''