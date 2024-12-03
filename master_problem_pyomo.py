# Import necessary libraries
from pyomo.environ import *
from pyomo.opt import SolverFactory
import pyomo.environ as pyo


# Initialize Benders Decomposition
def initial_master_problem(data, thermal_gens, renewable_gens, time_periods, gen_startup_categories, iteration):

    master = ConcreteModel()

    print('Building master problem')

    master.ug = Var(thermal_gens.keys(), time_periods.keys(), within=Binary) 
    master.vg = Var(thermal_gens.keys(), time_periods.keys(), within=Binary) 
    master.wg = Var(thermal_gens.keys(), time_periods.keys(), within=Binary) 
    master.dg = Var(((g,s,t) for g in thermal_gens for s in gen_startup_categories[g] for t in time_periods), within=Binary)

    #The beta variable that should be representing the sub problem cost
    master.beta = Var(bounds=(-1e9, None), within=Reals)

    #Master problem objective function
    master.obj = Objective(expr=sum(
                        sum(
                            thermal_gens[g]['piecewise_production'][0]['cost'] * master.ug[g, t]
                            + sum(thermal_gens[g]['startup'][s]['cost'] * master.dg[g, s, t] for s in gen_startup_categories[g])
                            for t in time_periods
                        )
                        for g in thermal_gens
                    ) + master.beta)

    #Constrains, only those with just binary variables, the constrains are 99% sure formulated correctly
    master.uptimet0 = Constraint(thermal_gens.keys())
    master.downtimet0 = Constraint(thermal_gens.keys())
    master.logicalt0 = Constraint(thermal_gens.keys())
    master.startupt0 = Constraint(thermal_gens.keys())
    master.shutdownt0 = Constraint(thermal_gens.keys())

    for g, gen in thermal_gens.items():
        if gen['unit_on_t0'] == 1:
            if iteration == -1:
                #print('first iteration-',g)
                for t in time_periods:
                #    master.ug[g,t] = 1  
                    master.uptimet0[g] = (master.ug[g,t] - 1) == 0 #(4)

            elif gen['time_up_minimum'] - gen['time_up_t0'] >= 1:
                master.uptimet0[g] = sum( (master.ug[g,t] - 1) for t in range(1, min(gen['time_up_minimum'] - gen['time_up_t0'], data['time_periods'])+1)) == 0 #(4)
        elif gen['unit_on_t0'] == 0:
            if gen['time_down_minimum'] - gen['time_down_t0'] >= 1:
                master.downtimet0[g] = sum( master.ug[g,t] for t in range(1, min(gen['time_down_minimum'] - gen['time_down_t0'], data['time_periods'])+1)) == 0 #(5)
        else:
            raise Exception('Invalid unit_on_t0 for generator {}, unit_on_t0={}'.format(g, gen['unit_on_t0']))

    master.logicalt0[g] = master.ug[g,1] - gen['unit_on_t0'] - master.vg[g,1] + master.wg[g,1] == 0 #(6)

    startup_expr = sum( 
                        sum( master.dg[g,s,t] 
                                for t in range(
                                                max(1, gen['startup'][s+1]['lag']-gen['time_down_t0']+1),
                                                min(gen['startup'][s+1]['lag']-1, data['time_periods'])+1
                                              )
                            )
                       for s,_ in enumerate(gen['startup'][:-1]))
    if isinstance(startup_expr, int):
        pass
    else:
        master.startupt0[g] = startup_expr == 0 #(7)

    shutdown_constr = gen['unit_on_t0']*(gen['power_output_t0']-gen['power_output_minimum']) <= gen['unit_on_t0']*(gen['power_output_maximum'] - gen['power_output_minimum']) - max((gen['power_output_maximum'] - gen['ramp_shutdown_limit']),0)*master.wg[g,1]

    if isinstance(shutdown_constr, bool):
        pass
    else:
        master.shutdownt0[g] = shutdown_constr  #(10)

    master.mustrun = Constraint(thermal_gens.keys(), time_periods.keys())
    master.logical = Constraint(thermal_gens.keys(), time_periods.keys())
    master.uptime = Constraint(thermal_gens.keys(), time_periods.keys())
    master.downtime = Constraint(thermal_gens.keys(), time_periods.keys())
    master.startup_select = Constraint(thermal_gens.keys(), time_periods.keys())

    for g, gen in thermal_gens.items():
        for t in time_periods:
            master.mustrun[g,t] = master.ug[g,t] >= gen['must_run'] #(11)

            if t > 1:
                master.logical[g,t] = master.ug[g,t] - master.ug[g,t-1] == master.vg[g,t] - master.wg[g,t] #(12)
            UT = min(gen['time_up_minimum'], data['time_periods'])
            if t >= UT:
                master.uptime[g,t] = sum(master.vg[g,t] for t in range(t-UT+1, t+1)) <= master.ug[g,t] #(13)
            DT = min(gen['time_down_minimum'], data['time_periods'])
            if t >= DT:
                master.downtime[g,t] = sum(master.wg[g,t] for t in range(t-DT+1, t+1)) <= 1-master.ug[g,t] #(14)
            master.startup_select[g,t] = master.vg[g,t] == sum(master.dg[g,s,t] for s,_ in enumerate(gen['startup'])) #(16)

    master.startup_allowed = Constraint(master.dg.index_set())
    for g, gen in thermal_gens.items():
        for s,_ in enumerate(gen['startup'][:-1]): ## all but last
            for t in time_periods:
                if t >= gen['startup'][s+1]['lag']:
                    master.startup_allowed[g,s,t] = master.dg[g,s,t] <= sum(master.wg[g,t-i] for i in range(gen['startup'][s]['lag'], gen['startup'][s+1]['lag'])) #(15)

    master.new = Constraint(time_periods.keys())
    for t,t_idx in time_periods.items():
        master.new[t] = sum( master.ug[g,t]*thermal_gens[g]['power_output_maximum'] for g in thermal_gens) + sum(renewable_gens[g]['power_output_maximum'][t_idx] for g in renewable_gens) >= data['demand'][t_idx] + data['reserves'][t_idx] #(NEW)


    #This is the way i add benders cuts, think it should be fine
    master.benders_cuts = ConstraintList()
    
    return master


def solving_master_problem(master, thermal_gens, renewable_gens, time_periods, gen_startup_categories, iteration):
    
    #solve the problem
    solver = SolverFactory('gurobi')
    solver.solve(master, options={'MIPGap': 0.01}, tee=False)

    master_var_values = {
            'ug': {(g, t): value(master.ug[g, t]) for g in thermal_gens.keys() for t in time_periods.keys()},
            'vg': {(g, t): value(master.vg[g, t]) for g in thermal_gens.keys() for t in time_periods.keys()},
            'wg': {(g, t): value(master.wg[g, t]) for g in thermal_gens.keys() for t in time_periods.keys()},
            'dg': {(g, s, t): value(master.dg[g, s, t]) for g in thermal_gens.keys() for s in gen_startup_categories[g] for t in time_periods.keys()}
        }

    master_obj_value = pyo.value(master.obj)
    beta = pyo.value(master.beta)

    return master_obj_value, master_var_values, beta
