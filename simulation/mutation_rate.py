#! /usr/local/bin/env python
# -*- coding: utf-8 -*-
import random
import numpy as np
import csv
import sys
import time
import copy

t1 = time.time()
# defined parameters
TSG_nonsyn_site = 191624
TSG_syn_site = 61470
control_nonsyn_site = 923307
control_syn_site = 319944
mean_age = 61.5
age_sd = 13.5


# define class of parameters
class parameter_object:
    def _init(self, N, generation_times, mutation_rate_coef, mutater_effect,
              s_exp_mean, selection_coef, TSGnon_s, mutater_s,
              syn_s=0, controlnon_s=0):
        self.N = N
        self.generation_times = generation_times
        self.mutation_rate = 1.5*(10**-8)*mutation_rate_coef
        self.mutater_effect = mutater_effect
        self.s_exp_mean = s_exp_mean
        self.selection_coef = selection_coef
        self.TSGnon_s = TSGnon_s
        self.mutater_s = mutater_s
        self.syn_s = syn_s
        self.controlnon_s = controlnon_s


def simulation(parameter_obj):
    tb = copy.copy(t1)
    cancer_population = Population(parameters)
    for t in range(parameters.generation_times):
        cancer_population.add_new_mutation()
        cancer_population.next_generation_wf(parameters)
        if t % (parameters.generation_times/100) == 0:
            print("now " + str(t) + " generation")
            t2 = time.time()
            elapsed_time = t2 - tb
            tb = t2
            print(f"spent: {elapsed_time}")

    # cancer_population.print_individuals(out_file)
    t2 = time.time()
    elapsed_time = t2 - t1
    print(f"all time spent: {elapsed_time}")


class Individual:
    def __init__(self, mutater, tsg_n, tsg_s, control_n, control_s, params):
        self._mutater = mutater
        self._tsg_n = tsg_n
        self._tsg_s = tsg_s
        self._control_n = control_n
        self._control_s = control_s
        age = mean_age
        age -=
        age += np.random.normal(0, age_sd)
        age = 100 if age > 100 else age
        age = 0 if age < 0 else age
        self._onset_age = age
        self._fitness = age / mean_age

    @property
    def fitness(self):
        return self._fitness

    @property
    def mutater(self):
        return self._mutater

    @property
    def onset_age(self):
        return self._onset_age

    # get [mutater, tsg_n, tsg_s, control_n, control_s]
    def mutations(self):
        mutations = copy.deepcopy([self._tsg_n, self._tsg_s])
        mutations += copy.deepcopy([self._control_n, self._control_s])
        return(mutations)

    def add_tsg_n(self, mutations):
        self._tsg_n = copy.copy(self._tsg_n) + (mutations)

    def add_tsg_s(self, mutations):
        self._tsg_s = copy.copy(self._tsg_s) + (mutations)

    def add_control_n(self, mutations):
        self._control_n = copy.copy(self._control_n) + (mutations)

    def add_control_s(self, mutations):
        self._control_s = copy.copy(self._control_s) + (mutations)

    def add_mutater(self, add_or_not):
        self._mutater = copy.copy(self.mutater) + add_or_not


# make de nove mutations list (list of each ind de novo mutation num)
def new_mutation(mp, site_num):
    new_mus = []
    for x in range(3):
        mutation_rate = 1.5*(10**(-8)) * (mutater_effect ** x) * site_num
        new_mus.extend(np.random.poisson(mutation_rate, mp[x]).tolist())
    return new_mus


# make offspring from two individuals
def reproduct(ind1, ind2):
    mutater = np.random.binomial(ind1.mutater + ind2.mutater, 0.5)
    mutater = 2 if mutater > 2 else mutater
    muts = [ind1.mutations()[i] + ind2.mutations()[i] for i in range(4)]
    new_mut = [random.sample(i, np.random.binomial(len(i), 0.5)) for i in muts]
    return(Individual(mutater, *new_mut))


class Population:
    def __init__(self, parameters):
        self.individuals = [Individual(mutater=0, tsg_n=[], tsg_s=[],
                                       control_n=[], control_s=[],
                                       params=parameters)
                            for i in range(parameters.N)]

    def get_fitness_list(self):
        fitness_list = [ind.fitness for ind in self.individuals]
        fit_sum = sum(fitness_list)
        fitness_list = [fit / fit_sum for fit in fitness_list]
        return fitness_list
    # add new mutations to each individuals

    def add_new_mutation(self, parameters):
        # individuals num of [mutater=0, mutater=1, mutater=2]
        muter_pnum = [len([x for x in self.individuals if x.mutater == 0])]
        muter_pnum.append(len([x for x in self.individuals if x.mutater == 1]))
        muter_pnum.append(len([x for x in self.individuals if x.mutater == 2]))
        mutater_mut_rate = parameters.mutation_rate*parameters.mutater_length
        new_mutater = np.random.binomial(1, mutater_mut_rate,
                                         parameters.N).tolist()
        new_mut_tn = new_mutation(muter_pnum, TSG_nonsyn_site)
        new_mut_ts = new_mutation(muter_pnum, TSG_syn_site)
        new_mut_cn = new_mutation(muter_pnum, control_nonsyn_site)
        new_mut_cs = new_mutation(muter_pnum, control_syn_site)
        for i in range(all_n):
            if self.individuals[i].mutater < 2:
                self.individuals[i].add_mutater(new_mutater[i])
            new_mus = np.random.randint(0, TSG_nonsyn_site,
                                        new_mut_tn[i]).tolist()
            self.individuals[i].add_tsg_n(new_mus)
            new_mus = np.random.randint(0, TSG_syn_site,
                                        new_mut_ts[i]).tolist()
            self.individuals[i].add_tsg_s(new_mus)
            new_mus = np.random.randint(0, control_nonsyn_site,
                                        new_mut_cn[i]).tolist()
            self.individuals[i].add_control_n(new_mus)
            new_mus = np.random.randint(0, control_syn_site,
                                        new_mut_cs[i]).tolist()
            self.individuals[i].add_control_s(new_mus)

    # make next generation population
    def next_generation_wf(self, population_N):
        fitness = self.get_fitness_list()
        rand_ind = np.random.choice(self.individuals, population_N*2,
                                    p=fitness)
        next_generation = [reproduct(rand_ind[i], rand_ind[i+1])
                           for i in range(0, population_N*2, 2)]
        self.individuals = next_generation
        self.individuals.sort(key=lambda ind: ind.mutater)

    def print_individuals(self, file):
        with open(file, "w") as f:
            out = csv.writer(f, delimiter="\t")
            col = ["id", "mutater", "TSG_nonsyn", "TSG_syn"]
            col += ["control_nonsyn", "control_syn", "onset_age"]
            col += ["fitness"]
            out.writerow(col)
            for row_n in range(len(self.individuals)):
                ind = self.individuals[row_n]
                mut = [len(muts) for muts in ind.mutations()]
                info = [ind.mutater, *mut]
                info += [ind.onset_age, ind.fitness]
                out.writerow([row_n, *info])
