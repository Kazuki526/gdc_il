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
tsg_non_site = 191624
tsg_syn_site = 61470
cont_non_site = 923307
cont_syn_site = 319944
mean_age = 61.5
age_sd = 13.5


# define class of parameters
class parameter_object:
    def _init(self, N, generation_times, mutation_rate_coef, mutater_effect,
              s_exp_mean, selection_coef, tsgnon_s, mutater_s,
              syn_s=0, contnon_s=0):
        self.N = N
        self.generation_times = generation_times
        self.mutation_rate = 1.5*(10**-8)*mutation_rate_coef
        self.mutater_effect = mutater_effect
        tsgnon_s = np.rondom.exponential(s_exp_mean, tsg_non_site).tolist()
        tsgnon_s = [1 if s > 1 else s for s in tsgnon_s]
        self.tsgnon_s = [tsgnon_s * round(exp, 4) for exp in tsgnon_s]
        tsgsyn_s = np.rondom.exponential(s_exp_mean, tsg_syn_site).tolist()
        tsgnon_s = [1 if s > 1 else s for s in tsgnon_s]
        self.tsgsyn_s = [syn_s * round(exp, 4) for exp in tsgsyn_s]
        contnon_s = np.rondom.exponential(s_exp_mean, cont_non_site).tolist()
        contnon_s = [1 if s > 1 else s for s in contnon_s]
        self.contnon_s = [contnon_s * round(exp, 4) for exp in contnon_s]
        contsyn_s = np.rondom.exponential(s_exp_mean, cont_syn_site).tolist()
        contnon_s = [1 if s > 1 else s for s in contnon_s]
        self.contsyn_s = [syn_s * round(exp, 4) for exp in contsyn_s]
        self.selection_coef = selection_coef
        self.mutater_s = mutater_s


def simulation(parameter_obj):
    tb = copy.copy(t1)
    cancer_population = Population(parameters=parameter_obj)
    for t in range(parameter_obj.generation_times):
        cancer_population.add_new_mutation()
        cancer_population.next_generation_wf(parameter_obj)
        if t % 100 == 0:
            t2 = time.time()
            elapsed_time = t2 - tb
            print(f"now {t} generation done, spent: {elapsed_time}")
            tb = t2
    # cancer_population.print_individuals(out_file)
    t2 = time.time()
    elapsed_time = t2 - t1
    print(f"all time spent: {elapsed_time}")


def genotype_divide(mutations):
    homo = [x for x in set(mutations) if mutations.count(x) > 1]
    hetero = list(set(mutations) - set(homo))
    return(homo, hetero)


class Individual:
    def __init__(self, mutater, tsg_non, tsg_syn, cont_non, cont_syn,
                 parameters):
        self._mutater = mutater
        self._tsg_non_hom, self._tsg_non_het = genotype_divide(tsg_non)
        self._tsg_syn_hom, self._tsg_syn_het = genotype_divide(tsg_syn)
        self._cont_non_hom, self.cont_non_het = genotype_divide(cont_non)
        self._cont_syn_hom, self.cont_syn_het = genotype_divide(cont_syn)
        age = mean_age
        age -= sum([parameters.tsgnon_s[mut] for mut in tsg_non])
        age -= sum([parameters.tsgsyn_s[mut] for mut in tsg_syn])
        age -= sum([parameters.contnon_s[mut] for mut in cont_non])
        age -= sum([parameters.contsyn_s[mut] for mut in cont_syn])
        age += np.random.normal(0, age_sd)
        age = 100 if age > 100 else age
        age = 20 if age < 20 else age
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

    # get [tsg_n_het, tsg_s_het, cont_n_het, cont_s_het]
    def mutations(self):
        mutations = copy.deepcopy([self._tsg_n_het, self._tsg_s_het])
        mutations += copy.deepcopy([self._cont_n_het, self._cont_s_het])
        return(mutations)

    # get [tsg_n_homo, tsg_s_homo, cont_n_homo, cont_n_homo]
    def homo_mutations(self):
        mutations = copy.deepcopy([self._tsg_n_hom, self._tsg_s_hom])
        mutations += copy.deepcopy([self._cont_n_hom, self._cont_s_hom])
        return(mutations)

    def add_tsg_non(self, muts):
        old_het = set(copy.copy(self._tsg_non_het))
        old_hom = set(copy.copy(self._tsg_non_hom))
        self._tsg_non_het = list(old_het - set(muts))+list(old_hom & set(muts))
        self._tsg_non_hom = list(old_hom - set(muts))+list(old_het & set(muts))

    def add_tsg_syn(self, muts):
        old_het = set(copy.copy(self._tsg_syn_het))
        old_hom = set(copy.copy(self._tsg_syn_hom))
        self._tsg_syn_het = list(old_het - set(muts))+list(old_hom & set(muts))
        self._tsg_syn_hom = list(old_hom - set(muts))+list(old_het & set(muts))

    def add_cont_non(self, muts):
        old_het = set(copy.copy(self._cont_n_het))
        old_hom = set(copy.copy(self._cont_n_hom))
        self._cont_non_het = list(old_het-set(muts))+list(old_hom & set(muts))
        self._cont_non_hom = list(old_hom-set(muts))+list(old_het & set(muts))

    def add_cont_syn(self, muts):
        old_het = set(copy.copy(self._cont_syn_het))
        old_hom = set(copy.copy(self._cont_syn_hom))
        self._cont_syn_het = list(old_het-set(muts))+list(old_hom & set(muts))
        self._cont_syn_hom = list(old_hom-set(muts))+list(old_het & set(muts))

    def add_mutater(self, add_or_not):
        self._mutater = copy.copy(self._mutater) + add_or_not


# make de nove mutations list (list of each ind de novo mutation num)
def new_mutation(mp, site_num, parameters):
    new_mus = []
    for x in range(3):
        mutation_rate = parameters.mutation_rate * site_num
        new_mus.extend(np.random.poisson(mutation_rate, mp[x]).tolist())
    return new_mus


# make offspring from two Individuals
def reproduct(ind1, ind2):
    mutater = np.random.binomial(ind1.mutater + ind2.mutater, 0.5)
    mutater = 2 if mutater > 2 else mutater
    muts = [ind1.mutations()[i] + ind2.mutations()[i] for i in range(4)]
    new_mut = [random.sample(l, np.random.binomial(len(l), 0.5)) for l in muts]
    new_mut = [new_mut[i] + ind1.homo_mutations[i] for i in range(4)]
    return(Individual(mutater, *new_mut))


class Population:
    def __init__(self, params):
        self.individuals = [Individual(mutater=0, tsg_non=[], tsg_syn=[],
                                       cont_non=[], control_syn=[],
                                       parametes=params)
                            for i in range(params.N)]

    def get_fitness_list(self):
        fitness_list = [ind.fitness for ind in self.individuals]
        fit_sum = sum(fitness_list)
        fitness_list = [fit / fit_sum for fit in fitness_list]
        return fitness_list

    # add new mutations to each individuals
    def add_new_mutation(self, params):
        # individuals num of [mutater=0, mutater=1, mutater=2]
        muter_pnum = [len([x for x in self.individuals if x.mutater == 0])]
        muter_pnum.append(len([x for x in self.individuals if x.mutater == 1]))
        muter_pnum.append(len([x for x in self.individuals if x.mutater == 2]))
        mutater_mut_rate = params.mutation_rate*params.mutater_length
        new_mutater = np.random.binomial(1, mutater_mut_rate,
                                         params.N).tolist()
        new_mut_tn = new_mutation(muter_pnum, tsg_non_site, params)
        new_mut_ts = new_mutation(muter_pnum, tsg_syn_site, params)
        new_mut_cn = new_mutation(muter_pnum, cont_non_site, params)
        new_mut_cs = new_mutation(muter_pnum, cont_syn_site, params)
        for n in range(params.N):
            if self.individuals[n].mutater < 2:
                self.individuals[n].add_mutater(new_mutater[n])
            if new_mut_tn[n] != 0:
                new_mus = np.random.randint(0, tsg_non_site,
                                            new_mut_tn[n]).tolist()
                self.individuals[n].add_tsg_n(new_mus)
            if new_mut_ts[n] != 0:
                new_mus = np.random.randint(0, tsg_syn_site,
                                            new_mut_ts[n]).tolist()
                self.individuals[n].add_tsg_s(new_mus)
            if new_mut_cn[n] != 0:
                new_mus = np.random.randint(0, cont_non_site,
                                            new_mut_cn[n]).tolist()
                self.individuals[n].add_cont_n(new_mus)
            if new_mut_cs[n] != 0:
                new_mus = np.random.randint(0, cont_syn_site,
                                            new_mut_cs[n]).tolist()
                self.individuals[n].add_cont_s(new_mus)

    # make next generation population
    def next_generation_wf(self, params):
        fitness = self.get_fitness_list()
        rand_ind = np.random.choice(self.individuals, params.N*2, p=fitness)
        next_generation = [reproduct(rand_ind[n], rand_ind[n+1])
                           for n in range(0, params.N*2, 2)]
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
