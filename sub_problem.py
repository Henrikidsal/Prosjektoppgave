# Import necessary libraries
from pyomo.environ import *
from pyomo.opt import SolverFactory
import pyomo.environ as pyo
from pyomo.environ import Suffix

# Subproblem
def subproblem(data, master_var_values, thermal_gens, renewable_gens, time_periods, gen_pwl_points):

    m = ConcreteModel()

    #this is for collecting dual variables
    m.dual = Suffix(direction=Suffix.IMPORT)

    #sub problem variables
    m.cg = Var(thermal_gens.keys(), time_periods.keys())
    m.pg = Var(thermal_gens.keys(), time_periods.keys(), within=NonNegativeReals)  
    m.rg = Var(thermal_gens.keys(), time_periods.keys(), within=NonNegativeReals)  
    m.pw = Var(renewable_gens.keys(), time_periods.keys(), within=NonNegativeReals)
    m.lg = Var(((g,l,t) for g in thermal_gens for l in gen_pwl_points[g] for t in time_periods), within=UnitInterval) ##
    
    #Slack variables to make the problem feasable
    m.slack_demand = Var(time_periods.keys(), within=NonNegativeReals)
    #m.slack_reserve = Var(time_periods.keys(), within=NonNegativeReals)
    PENALTY = 400

    #The variables from the master problem, an alternative way of using Variables that i fix is below
    m.ug = Var(thermal_gens.keys(), time_periods.keys(), within=NonNegativeReals)
    m.vg = Var(thermal_gens.keys(), time_periods.keys(), within=NonNegativeReals)
    m.wg = Var(thermal_gens.keys(), time_periods.keys(), within=NonNegativeReals)

    '''
    # Fix variables to master solution values
    for g in thermal_gens.keys():
        for t in time_periods.keys():
            m.ug[g, t].fix(master_solution['ug'][g, t])
            m.vg[g, t].fix(master_solution['vg'][g, t])
            m.wg[g, t].fix(master_solution['wg'][g, t])
    '''
    # Objective Function 
    m.obj = Objective(expr= sum(sum(m.cg[g, t] for t in time_periods) for g in thermal_gens) + PENALTY * (sum(m.slack_demand[t] for t in time_periods) ))

    # Constraints (Pretty sure these are correct, since they are the same as in a different formulation without benders)
    m.fixing_ug = Constraint(thermal_gens.keys(), time_periods.keys())
    m.fixing_vg = Constraint(thermal_gens.keys(), time_periods.keys())
    m.fixing_wg = Constraint(thermal_gens.keys(), time_periods.keys())

    for g in thermal_gens.keys():
        for t in time_periods.keys():
            m.fixing_ug[g, t] = m.ug[g, t] == master_var_values['ug'][g, t]
            m.fixing_vg[g, t] = m.vg[g, t] == master_var_values['vg'][g, t]
            m.fixing_wg[g, t] = m.wg[g, t] == master_var_values['wg'][g, t]

    m.demand = Constraint(time_periods.keys())
    m.reserves = Constraint(time_periods.keys())
    for t,t_idx in time_periods.items():
        m.demand[t] = sum( m.pg[g,t]+gen['power_output_minimum']*m.ug[g,t] for (g, gen) in thermal_gens.items() ) + sum(m.pw[w,t] for w in renewable_gens) + m.slack_demand[t] == data['demand'][t_idx] #(2)
        m.reserves[t] = sum( m.rg[g,t] for g in thermal_gens ) >= data['reserves'][t_idx] #(3)

    m.rampupt0 = Constraint(thermal_gens.keys())
    m.rampdownt0 = Constraint(thermal_gens.keys())
    for g, gen in thermal_gens.items():
        m.rampupt0[g] = m.pg[g,1] + m.rg[g,1] - gen['unit_on_t0']*(gen['power_output_t0']-gen['power_output_minimum']) <= gen['ramp_up_limit'] #(8)
        m.rampdownt0[g] = gen['unit_on_t0']*(gen['power_output_t0']-gen['power_output_minimum']) - m.pg[g,1] <= gen['ramp_down_limit'] #(9)

    m.gen_limit1 = Constraint(thermal_gens.keys(), time_periods.keys())
    m.gen_limit2 = Constraint(thermal_gens.keys(), time_periods.keys())
    m.ramp_up = Constraint(thermal_gens.keys(), time_periods.keys())
    m.ramp_down = Constraint(thermal_gens.keys(), time_periods.keys())
    m.power_select = Constraint(thermal_gens.keys(), time_periods.keys())
    m.cost_select = Constraint(thermal_gens.keys(), time_periods.keys())
    m.on_select = Constraint(thermal_gens.keys(), time_periods.keys())

    for g, gen in thermal_gens.items():
        for t in time_periods:
            m.gen_limit1[g,t] = m.pg[g,t]+m.rg[g,t]+ max((gen['power_output_maximum'] - gen['ramp_startup_limit']),0)*m.vg[g,t] - (gen['power_output_maximum'] - gen['power_output_minimum'])*m.ug[g,t] <= 0 #(17)
            if t < len(time_periods): 
                m.gen_limit2[g,t] = m.pg[g,t]+m.rg[g,t] + max((gen['power_output_maximum'] - gen['ramp_shutdown_limit']),0)*m.wg[g,t+1] - (gen['power_output_maximum'] - gen['power_output_minimum'])*m.ug[g,t] <= 0 #(18)
            if t > 1:
                m.ramp_up[g,t] = m.pg[g,t]+m.rg[g,t] - m.pg[g,t-1] <= gen['ramp_up_limit'] #(19)
                m.ramp_down[g,t] = m.pg[g,t-1] - m.pg[g,t] <= gen['ramp_down_limit'] #(20)
            piece_mw1 = gen['piecewise_production'][0]['mw']
            piece_cost1 = gen['piecewise_production'][0]['cost']
            m.power_select[g,t] = m.pg[g,t] == sum( (piece['mw'] - piece_mw1)*m.lg[g,l,t] for l,piece in enumerate(gen['piecewise_production'])) #(21)
            m.cost_select[g,t] = m.cg[g,t] == sum( (piece['cost'] - piece_cost1)*m.lg[g,l,t] for l,piece in enumerate(gen['piecewise_production'])) #(22)
            m.on_select[g,t] = m.ug[g,t] == sum(m.lg[g,l,t] for l,_ in enumerate(gen['piecewise_production'])) #(23)

    for w, gen in renewable_gens.items():
        for t, t_idx in time_periods.items():
            m.pw[w,t].setlb(gen['power_output_minimum'][t_idx]) #(24)
            m.pw[w,t].setub(gen['power_output_maximum'][t_idx]) #(24)

    #solves sub problem
    solver = SolverFactory('gurobi')
    #solver.options['OutputFlag'] = 0
    #solver.options['DualReductions'] = 0
    solver.solve(m, tee=False)
    theta_j = pyo.value(m.obj)

    total_slack = sum(pyo.value(m.slack_demand[t]) for t in time_periods)
    if total_slack > 1e-6:  # Use a small threshold to account for numerical precision
        print("infeasable")
    else:
        print("FEASABLE")

    
    dual_values = {
        "fixing_ug": {},
        "fixing_vg": {},
        "fixing_wg": {}
    }

    # Loop over thermal generators and time periods to populate the duals
    for g in thermal_gens.keys():
        for t in time_periods.keys():
            # Save the dual value for fixing_ug
            dual_values["fixing_ug"][(g, t)] = m.dual[m.fixing_ug[g, t]]
            
            # Save the dual value for fixing_vg
            dual_values["fixing_vg"][(g, t)] = m.dual[m.fixing_vg[g, t]]
            
            # Save the dual value for fixing_wg
            dual_values["fixing_wg"][(g, t)] = m.dual[m.fixing_wg[g, t]]
    
    # Print the values of the slack variables
    print("\nSlack Variable Values:") 
    sum_slack_demand = 0
    for t in time_periods.keys():
        slack_demand_value = pyo.value(m.slack_demand[t])
        sum_slack_demand += slack_demand_value
        #slack_reserve_value = pyo.value(m.slack_reserve[t])
        #print(f"Time Period {t}: Slack Demand = {slack_demand_value}")
    print(f"total slack = to {sum_slack_demand}")

    #dual values goes to cuts, theta_j goes to master for creating UB
    return dual_values, theta_j


'''
#This is collecting the dual values in a dictionary
dual_values = {}
relevant_constraints = [m.fixing_ug, m.fixing_vg, m.fixing_wg]
try:
    for constraint in relevant_constraints:
        constr_name = constraint.local_name
        for index in constraint:
            con = constraint[index]
            dual = m.dual.get(con, 0)
            dual_values[(constr_name, index)] = dual
except KeyError as e:
    print(f"Error accessing duals: {e}")
'''