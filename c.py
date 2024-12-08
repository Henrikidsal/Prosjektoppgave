from pyomo.environ import *
import json
import sys
import timeit
import numpy as np
import os.path
from pathlib import Path
import itertools
import dimod
import sqlite3
import tabu
import pprint
from dwave.system import LeapHybridCQMSampler
from dimod import BinaryQuadraticModel
from tqdm import tqdm
from dimod import quicksum as Sum
import dimod
import time
import pandas as pd

## Grab instance file from first command line argument
data_file = sys.argv[1]

print('loading data')
data = json.load(open(data_file, 'r'))
print(data['demand'][1])

thermal_gens = data['thermal_generators']
renewable_gens = data['renewable_generators']


time_periods = {t+1 : t for t in range(data['time_periods'])}


gen_startup_categories = {g : list(range(0, len(gen['startup']))) for (g, gen) in thermal_gens.items()}
#gen_pwl_points = {g : list(range(0, 1)) for (g, gen) in thermal_gens.items()}
gen_pwl_points = {g : list(range(0, len(gen['piecewise_production']))) for (g, gen) in thermal_gens.items()}
#gen_pwl_points = gen_startup_categories
#gen_pwl_points = gen_startup_categories
print('building model')
#m = ConcreteModel()
m = dimod.ConstrainedQuadraticModel()

#for (g, gen) in thermal_gens.items():
    #print(gen['piecewise_production'])
    #print(gen['piecewise_production'][0]['mw'])
    #print(gen_pwl_points[g])#

#print(gen_pwl_points)

# Assuming you have already defined your QUBO model
#cqm = dimod.ConstrainedQuadraticModel()

pg, rg, cg, pw, ug, vg, wg, dg, lg = {}, {}, {}, {}, {}, {}, {}, {}, {}

#print(time_periods.keys())
time_periods = dict(itertools.islice(time_periods.items(), 2))

#thermal_gens = dict(itertools.islice(thermal_gens.items(), 4))
#gen_pwl_points = dict(itertools.islice(gen_pwl_points.items(), 4))
#print(gen_pwl_points.keys())
#print(thermal_gens['115_STEAM_1'])
#print(gen_pwl_points['115_STEAM_1'])


for g in thermal_gens.keys():
    #g_0 = 
    for t in time_periods.keys():
                                    
        cg[g, t] = dimod.Real(('cg',g,t),lower_bound= 0)
        pg[g, t] = dimod.Real(('pg',g,t),lower_bound = 0)
        rg[g, t] = dimod.Real(('rg',g,t), lower_bound = 0)
        
        
        ug[g, t] = dimod.Binary(('ug',g,t))
        vg[g, t] = dimod.Binary(('vg',g,t))
        wg[g, t] = dimod.Binary(('wg',g,t))
        
#print(renewable_gens.keys())      

for w, gen in renewable_gens.items():
    for t, t_idx in time_periods.items():
        pw[w, t] = dimod.Real(('pw',w,t),lower_bound=gen['power_output_minimum'][t_idx],upper_bound=gen['power_output_maximum'][t_idx])
"""
for w in renewable_gens.keys():
    for t  in time_periods.keys():
        pw[w, t] = dimod.Real(('pw',w,t),lower_bound=0)

"""

for g in thermal_gens:
    for s in gen_startup_categories[g]:
        for t in time_periods:
            dg[g,s,t] = dimod.Binary(('dg',g,s,t))
    for l in range(1):#gen_pwl_points[g]:
        #print(l)
        for t in time_periods:
            lg[g,l,t] = dimod.Real(('lg',g,l,t),lower_bound = 0, upper_bound=1)

print("Here",len(lg))  



m.set_objective(sum(
                          sum(cg[g,t] + gen['piecewise_production'][0]['cost']*ug[g,t]
                              + sum( gen_startup['cost']*dg[g,s,t] for (s, gen_startup) in enumerate(gen['startup']))
                          for t in time_periods)
                        for g, gen in thermal_gens.items() )
                        )


for t,t_idx in time_periods.items():
    m.add_constraint(sum(pg[g,t] + gen['power_output_minimum']*ug[g,t] for (g, gen) in thermal_gens.items() ) + sum(pw[w,t] for w in renewable_gens ) == data['demand'][t_idx], "demand_%s" % t)#(2)-------------------_> WHY THIS
    m.add_constraint(sum( rg[g,t] for g in thermal_gens ) >= data['reserves'][t_idx], "reserves_%s" % t) #(3)



for g, gen in thermal_gens.items():
    if gen['unit_on_t0'] == 1:
        if gen['time_up_minimum'] - gen['time_up_t0'] >= 1:
            m.add_constraint(sum( (ug[g,t] - 1) for t in range(1, min(gen['time_up_minimum'] - gen['time_up_t0'], data['time_periods'])+1)) == 0, "uptimet0_%s" % g) #(4)
            
            
    elif gen['unit_on_t0'] == 0:
        if gen['time_down_minimum'] - gen['time_down_t0'] >= 1:
            m.add_constraint(sum( ug[g,t] for t in range(1, min(gen['time_down_minimum'] - gen['time_down_t0'], data['time_periods'])+1)) == 0, "downtimet0_%s" % g )#(5)
    else:
        raise Exception('Invalid unit_on_t0 for generator {}, unit_on_t0={}'.format(g, gen['unit_on_t0']))

    m.add_constraint(ug[g,1] - gen['unit_on_t0'] - vg[g,1] + wg[g,1] == 0, "logicalt0%s" % g) #(6)


    startup_expr = sum( 
                        sum( dg[g,s,t] 
                                for t in range(
                                                max(1,gen['startup'][s+1]['lag']-gen['time_down_t0']+1),
                                                min(gen['startup'][s+1]['lag']-1,len(time_periods))+1
                                              )
                            ) 
                       for s,_ in enumerate(gen['startup'][:-1])) ## all but last
    if isinstance(startup_expr, int):
        pass
    else:
        m.add_constraint(sum( 
                        sum( dg[g,s,t] 
                                for t in range(
                                                max(1,gen['startup'][s+1]['lag']-gen['time_down_t0']+1),
                                                min(gen['startup'][s+1]['lag']-1,len(time_periods))+1
                                              )
                            ) 
                       for s,_ in enumerate(gen['startup'][:-1])) == 0, "startupt0_%s_%s_%s" % (g,s,t))#(7)

    m.add_constraint(pg[g,1] + rg[g,1] - gen['unit_on_t0']*(gen['power_output_t0']-gen['power_output_minimum']) <= gen['ramp_up_limit'],"rampupt0_%s" % g )#(8)

    m.add_constraint(gen['unit_on_t0']*(gen['power_output_t0']-gen['power_output_minimum']) - pg[g,1] <= gen['ramp_down_limit'], "rampdownt0_%s" % g)#(9)


    shutdown_constr = gen['unit_on_t0']*(gen['power_output_t0']-gen['power_output_minimum']) <= gen['unit_on_t0']*(gen['power_output_maximum'] - gen['power_output_minimum']) - max((gen['power_output_maximum'] - gen['ramp_shutdown_limit']),0)*wg[g,1] #(10)

    if isinstance(shutdown_constr, bool):
        pass
    else:
        m.add_constraint(gen['unit_on_t0']*(gen['power_output_t0']-gen['power_output_minimum']) <= gen['unit_on_t0']*(gen['power_output_maximum'] - gen['power_output_minimum']) - max((gen['power_output_maximum'] - gen['ramp_shutdown_limit']),0)*wg[g,1],"shutdownt0_%s" % g)








for g, gen in thermal_gens.items():
    for t in time_periods:
        m.add_constraint(ug[g,t] >= gen['must_run'],"mustrun_%s_%s" % (g, t)) #(11)

        if t > 1:
            m.add_constraint(ug[g,t] - ug[g,t-1] - vg[g,t] + wg[g,t] == 0,"logical_%s_%s" % (g,t)) #(12)

        UT = min(gen['time_up_minimum'],data['time_periods'])
        if t >= UT:
            m.add_constraint(sum(vg[g,t] for t in range(t-UT+1, t+1)) - ug[g,t] <= 0, "uptime_%s_%s" % (g,t)) #(13)

        DT = min(gen['time_down_minimum'],data['time_periods'])
        if t >= DT:
            m.add_constraint(sum(wg[g,t] for t in range(t-DT+1, t+1))  + ug[g,t] <= 1, "downtime_%s_%s" % (g,t))#(14)
        m.add_constraint(vg[g,t] - sum(dg[g,s,t] for s,_ in enumerate(gen['startup'])) == 0, "startup_select_%s_%s" % (g,t) ) #(16)

        m.add_constraint(pg[g,t] + rg[g,t] - (gen['power_output_maximum'] - gen['power_output_minimum']) * ug[g,t] + max((gen['power_output_maximum'] - gen['ramp_startup_limit']),0) * vg[g,t] <= 0, "gen_limit1_%s_%s" % (g,t))#(17)-------> and also this one

        if t < len(time_periods): 
            m.add_constraint(pg[g,t] + rg[g,t] - (gen['power_output_maximum'] - gen['power_output_minimum']) * ug[g,t] + max((gen['power_output_maximum'] - gen['ramp_shutdown_limit']),0) * wg[g,t+1] <= 0, "gen_limit2_%s_%s" % (g,t))#(18)#---------__> this one!!!

        if t > 1:
            m.add_constraint(pg[g,t] + rg[g,t] - pg[g,t-1] <=  gen['ramp_up_limit'], "ramp_up_%s_%s" % (g,t)) #(19)
            m.add_constraint(pg[g,t-1] - pg[g,t] <= gen['ramp_down_limit'], "ramp_down_%s_%s" % (g,t)) #(20

        piece_mw1 = gen['piecewise_production'][0]['mw']
        piece_cost1 = gen['piecewise_production'][0]['cost']
        """ 
        m.add_constraint(pg[g,t] - sum( (piece['mw'] - piece_mw1) * lg[g,1,t] for l,piece in enumerate(gen['piecewise_production'])) == 0, "power_select_%s_%s" % (g,t)) #(21)
        m.add_constraint(cg[g,t] - sum( (piece['cost'] - piece_cost1) * lg[g,l,t] for l,piece in enumerate(gen['piecewise_production'])) == 0, "cost_select_%s_%s" % (g,t)) #(22)
        m.add_constraint(ug[g,t] - sum(lg[g,l,t] for l,_ in enumerate(gen['piecewise_production'])) == 0, "on_select_%s_%s" % (g,t)) #(23)
        """
                
        m.add_constraint(pg[g,t] - ( (gen['piecewise_production'][0]['mw'] - piece_mw1)) * lg[g,0,t]  == 0, "power_select_%s_%s" % (g,t)) #(21)
        m.add_constraint(cg[g,t] - ( (gen['piecewise_production'][0]['cost'] - piece_cost1)) * lg[g,0,t] == 0, "cost_select_%s_%s" % (g,t)) #(22)
        m.add_constraint(ug[g,t] - (lg[g,0,t] ) == 0, "on_select_%s_%s" % (g,t)) #(23)

for g, gen in thermal_gens.items():
    for s,_ in enumerate(gen['startup'][:-1]): ## all but last
        for t in time_periods:
            if t >= gen['startup'][s+1]['lag']:
                m.add_constraint(dg[g,s,t] - sum(wg[g,t-i] for i in range(gen['startup'][s]['lag'], gen['startup'][s+1]['lag'])) <= 0,"startup_allowed_%s_%s_%s" % (g,s,t) )#(15)

print("Total number of constraints:", len(m.constraints))
"""

for w, gen in renewable_gens.items():
    for t, t_idx in time_periods.items():
        print(w)
        m.set_lower_bound(('pw', w, t), gen['power_output_minimum'][t_idx])
        m.set_upper_bound(('pw', w, t), gen['power_output_maximum'][t_idx])
        #m.pw[w,t].setlb(gen['power_output_minimum'][t_idx]) #(24)
        #m.pw[w,t].setub(gen['power_output_maximum'][t_idx]) #(24)

"""
print("model setup complete")
def parse_best(sampleset):#find best solution from feasible sampleset
    print("parsing best ...")
    best = sampleset.filter(lambda row: row.is_feasible).first
    #s = [val for key, val in best.sample.items() if "s_" in key]
    #p = [val for key, val in best.sample.items() if "p_" in key]
    #r = [p*s for p, s in zip(p, s)]
    print("found best")
    return best



num_variables = len(m.variables)
num_constraints = len(m.constraints)

#m.check_feasible()
count = 0
for i in m.constraints:
    if 'startupt0' in i:
        count +=1
        #print(i)
        #pprint(i)
    else:
        continue
    #print(i)

print("Constraint count:",count)



print(f"Total number of variables: {num_variables}")
print(f"Total number of constraints: {num_constraints}")

#print("Total number of variables",len(cg)+len(pg)+len(rg)+len(lg)+len(vg)+len(wg)+len(dg)+len(pw)+len(ug))
#exit()
sampler = LeapHybridCQMSampler(profile='bocap')
#sampler.properties['minimum_it_s'] = 10
#print(sampler.properties['minimum_time_limit_s'])

#sampler = dimod.ExactCQMSolver()
#sampler = tabu.TabuSampler()
runs = 70
for i in range(runs):
    print("Run number:",i+1)

    if i == 0:
        start = timeit.default_timer()
        first_sampleset = sampler.sample_cqm(m,label="Unit commitment")
        #sampleset = sampler.sample(m)
        stop = timeit.default_timer()
        #first_sampleset = first_sampleset.to_serializable()
        run_time = first_sampleset.info['run_time']
        print("RUN TIME",run_time)
       
    elif i == 1:
        start = timeit.default_timer()
        sampleset = sampler.sample_cqm(m,label="Unit commitment")
        #sampleset = sampler.sample(m)
        stop = timeit.default_timer()
        #sampleset = sampleset.to_serializable()
        larger_sampleset = dimod.concatenate((first_sampleset,sampleset))
    else:
        start = timeit.default_timer()
        sampleset = sampler.sample_cqm(m,label="Unit commitment")
        #sampleset = sampler.sample(m)
        stop = timeit.default_timer()
        #sampleset = sampleset.to_serializable()
        larger_sampleset = dimod.concatenate((larger_sampleset,sampleset))
        
    
        




general = stop - start
print("Solve time:",general)


"""
with open('/home/quintonf/miniconda3/envs/CondaDWave/unit_commitment/solutions/data_dwave.json', 'w') as f:
    # write the dictionary to the file in JSON format
    json.dump(sampleset_series, f)

"""




#sampleset_pd = pd.DataFrame(sampleset)
#sampleset = dimod.drop_variables(sampleset)
#print(sampleset.first)
#print("info",sampleset.info)
#print("record",sampleset.record)
#print('record energy',sampleset.record.energy)
#print("Is feasible?",m.check_feasible(larger_sampleset))
sampleset_series = larger_sampleset.to_serializable()

#sampleset_series = first_sampleset.to_serializable()
#print(sampleset_series.record)
#sampleset_rec = sampleset.record
#sampleset_pd = pd.DataFrame(sampleset_rec)

with open('/home/quintonf/miniconda3/envs/CondaDWave/unit_commitment/solutions/data_dwave_2h_reduced_70_runs_final.json', 'w') as f:
    # write the dictionary to the file in JSON format
    json.dump(sampleset_series, f)
print("Saved sampleset")
best_one = parse_best(larger_sampleset)



print(best_one.energy)

exit()
sampleset_energy_pd.to_csv('/home/quintonf/miniconda3/envs/CondaDWave/unit_commitment/solutions/dwave_energy_solution.csv')


sampleset_energy_pd = sampleset_energy.to_pandas_dataframe()
print(sampleset_energy_pd)
solution = best_one.energy
#np.savetxt(solution,'/home/quintonf/miniconda3/envs/CondaDWave/unit_commitment/solutions/gurobi_rts_gmlc.csv',delimiter=',')
exit()
from pyomo.opt import SolverFactory
cbc = SolverFactory('gurobi')

print("solving")
t1=timeit.default_timer()
cbc.solve(m, options={'ratioGap':0.01}, tee=True)
time = timeit.default_timer()-t1

solution = pyomo.environ.value(m.obj)
s_path = str("/home/quintonf/miniconda3/envs/CondaDWave/unit_commitment/solutions/")
np.savetxt(solution,os.path.join(s_path,sol_file),delimiter=',')