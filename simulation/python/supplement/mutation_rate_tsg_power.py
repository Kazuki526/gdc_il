#! /usr/local/bin/env python
# -*- coding: utf-8 -*-
import random
import numpy as np
import time
import copy
from sklearn import linear_model
import csv

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
    def __init__(self, N, mutation_rate_coef, mutater_effect, s_exp_mean,
                 selection_coef, tsgnon_s, mutater_s, syn_s=0,
                 mutater_length=10000):
        self.N = N
        self.mutation_rate = 1.5*(10**-8)*mutation_rate_coef
        self.mutater_effect = mutater_effect
        self.mutater_length = mutater_length
        tsgnon_s_l = np.random.exponential(s_exp_mean, tsg_non_site).tolist()
        tsgnon_s_l = [1 if s > 1 else s for s in tsgnon_s_l]
        self.tsgnon_s = [round(tsgnon_s * s, 4) for s in tsgnon_s_l]
        tsgsyn_s_l = np.random.exponential(s_exp_mean, tsg_syn_site).tolist()
        tsgsyn_s_l = [1 if s > 1 else s for s in tsgsyn_s_l]
        self.tsgsyn_s = [round(syn_s * s, 4) for s in tsgsyn_s_l]
        self.selection_coef = selection_coef
        self.mutater_s = mutater_s


def genotype_divide(mutations):
    homo = [x for x in set(mutations) if mutations.count(x) > 1]
    hetero = list(set(mutations) - set(homo))
    return(homo, hetero)


class Individual:
    def __init__(self, mutater, tsg_non, tsg_syn, parameters):
        self._mutater = mutater
        self._tsg_non_hom, self._tsg_non_het = genotype_divide(tsg_non)
        self._tsg_syn_hom, self._tsg_syn_het = genotype_divide(tsg_syn)
        age = mean_age
        age -= sum([parameters.tsgnon_s[mut] for mut in tsg_non])
        age -= sum([parameters.tsgsyn_s[mut] for mut in tsg_syn])
        age += np.random.normal(0, age_sd)
        age = 100 if age > 100 else age
        age = 0 if age < 0 else age
        self._onset_age = age
        self._fitness = 1 - (1 - age / mean_age) * parameters.selection_coef

    @property
    def fitness(self):
        return self._fitness

    @property
    def mutater(self):
        return self._mutater

    @property
    def onset_age(self):
        return self._onset_age

    # get [tsg_non_het, tsg_syn_het]
    def mutations(self):
        mutations = copy.deepcopy([self._tsg_non_het, self._tsg_syn_het])
        return(mutations)

    # get [tsg_non_homo, tsg_syn_homo]
    def homo_mutations(self):
        mutations = copy.deepcopy([self._tsg_non_hom, self._tsg_syn_hom])
        return(mutations)

    def add_tsg_non(self, muts):
        old_het = set(copy.copy(self._tsg_non_het))
        old_hom = set(copy.copy(self._tsg_non_hom))
        self._tsg_non_het = list(old_het ^ set(muts))+list(old_hom & set(muts))
        self._tsg_non_hom = list(old_hom ^ set(muts))+list(old_het & set(muts))

    def add_tsg_syn(self, muts):
        old_het = set(copy.copy(self._tsg_syn_het))
        old_hom = set(copy.copy(self._tsg_syn_hom))
        self._tsg_syn_het = list(old_het ^ set(muts))+list(old_hom & set(muts))
        self._tsg_syn_hom = list(old_hom ^ set(muts))+list(old_het & set(muts))

    def add_mutater(self, add_or_not):
        self._mutater = copy.copy(self._mutater) + add_or_not

    def variant_num_tsg_non(self, variant_list):
        num = len(set(self._tsg_non_het) & set(variant_list))
        num += len(set(self._tsg_non_hom) & set(variant_list)) * 2
        return(num)

    def variant_num_tsg_syn(self, variant_list):
        num = len(set(self._tsg_syn_het) & set(variant_list))
        num += len(set(self._tsg_syn_hom) & set(variant_list)) * 2
        return(num)


# make de nove mutations list (list of each ind de novo mutation num)
def new_mutation(mp, site_num, params):
    new_mus = []
    for x in range(3):
        mutation_rate = params.mutation_rate * site_num
        new_mus.extend(np.random.poisson(mutation_rate, mp[x]).tolist())
    return new_mus


# make offspring from two Individuals
def reproduct(ind1, ind2, params):
    mutater = np.random.binomial(ind1.mutater + ind2.mutater, 0.5)
    mutater = 2 if mutater > 2 else mutater
    muts = [ind1.mutations()[i] + ind2.mutations()[i] for i in range(2)]
    new_mut = [random.sample(l, np.random.binomial(len(l), 0.5)) for l in muts]
    new_mut = [new_mut[i] + ind1.homo_mutations()[i] for i in range(2)]
    return(Individual(mutater, *new_mut, params))


class Population:
    def __init__(self, params):
        self.individuals = [Individual(mutater=0, tsg_non=[], tsg_syn=[],
                                       parameters=params)
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
        for n in range(params.N):
            if self.individuals[n].mutater < 2:
                self.individuals[n].add_mutater(new_mutater[n])
            if new_mut_tn[n] != 0:
                new_mus = np.random.randint(0, tsg_non_site,
                                            new_mut_tn[n]).tolist()
                self.individuals[n].add_tsg_non(new_mus)
            if new_mut_ts[n] != 0:
                new_mus = np.random.randint(0, tsg_syn_site,
                                            new_mut_ts[n]).tolist()
                self.individuals[n].add_tsg_syn(new_mus)

    # make next generation population
    def next_generation_wf(self, params):
        fitness = self.get_fitness_list()
        rand_ind = np.random.choice(self.individuals, params.N*2, p=fitness)
        next_generation = [reproduct(rand_ind[n], rand_ind[n+1], params)
                           for n in range(0, params.N*2, 2)]
        self.individuals = next_generation
        self.individuals.sort(key=lambda ind: ind.mutater)

    def variant_allelecount(self):
        v_AC = [[0 for i in range(tsg_non_site)]]
        v_AC += [[0 for i in range(tsg_syn_site)]]
        for ind in self.individuals:
            het_muts = ind.mutations()
            hom_muts = ind.homo_mutations()
            for i in range(2):
                for het_mut in het_muts[i]:
                    v_AC[i][het_mut] += 1
                for hom_mut in hom_muts[i]:
                    v_AC[i][hom_mut] += 2
        return(v_AC)

    def print_rare_variant_num(self, params):
        v_AC = self.variant_allelecount()
        rare_nums = [0, 0]
        for i in range(2):
            rare_nums[i] = sum([v_num for v_num in v_AC[i]
                                if v_num < params.N*2*0.0005])
        return(rare_nums)

    def print_summary(self, params, sample_num=6000):
        sample_num = params.N if sample_num > params.N else sample_num
        v_AC = self.variant_allelecount()
        rare_variants = [[], []]
        rare_nums = []
        for i in range(2):
            rare_variants[i] = [v for v in range(len(v_AC[i]))
                                if v_AC[i][v] < params.N*2*0.0005]
            rare_nums += [sum([v_AC[i][v] for v in rare_variants[i]])/params.N]
        sample = random.sample(self.individuals, sample_num)
        sample_age = [ind.onset_age for ind in sample]
        sample_tn_num = [ind.variant_num_tsg_non(rare_variants[0])
                         for ind in sample]
        sample_ts_num = [ind.variant_num_tsg_syn(rare_variants[1])
                         for ind in sample]
        sample_age = np.array(sample_age).reshape(-1, 1)
        sample_tn_num = np.array(sample_tn_num).reshape(-1, 1)
        sample_ts_num = np.array(sample_ts_num).reshape(-1, 1)
        tsg_non_reg = linear_model.LinearRegression()
        tsg_non_reg.fit(sample_age, sample_tn_num)
        tsg_syn_reg = linear_model.LinearRegression()
        tsg_syn_reg.fit(sample_age, sample_ts_num)
        return(*rare_nums, float(tsg_non_reg.coef_), float(tsg_syn_reg.coef_))

    def print_individuals(self, params, file="simulation_result.tsv"):
        v_AC = self.variant_allelecount()
        with open(file, "w") as f:
            out = csv.writer(f, delimiter="\t")
            col = ["id", "mutater", "TSG_nonsyn_het", "TSG_syn_het"]
            col += ["TSG_nonsyn_hom", "TSG_syn_hom", "onset_age"]
            col += ["fitness"]
            out.writerow(col)
            for row_n in range(len(self.individuals)):
                ind = self.individuals[row_n]
                out_mut_het, out_mut_hom = [[], []], [[], []]
                mut_het, mut_hom = ind.mutations(), ind.homo_mutations()
                for i in range(2):
                    out_mut_het[i] = [v for v in mut_het[i]
                                      if v_AC[i][v] < params.N*2*0.0005]
                    out_mut_hom[i] = [v for v in mut_hom[i]
                                      if v_AC[i][v] < params.N*2*0.0005]
                info = [ind.mutater, *out_mut_het, *out_mut_hom]
                info += [ind.onset_age, ind.fitness]
                out.writerow([row_n, *info])


def simulation(parameter_obj):
    population = Population(params=parameter_obj)
    focal = True
    v_num_lists = [[], []]
    generation = 0
    tb = copy.copy(t1)
    while focal:
        generation += 1
        population.add_new_mutation(parameter_obj)
        population.next_generation_wf(parameter_obj)
        if generation % 10 == 0:
            v_num_lists_sorted = [sorted(v_num_lists[i]) for i in range(2)]
            v_nums = population.print_rare_variant_num(parameter_obj)
            if len(v_num_lists[0]) < 10:
                for i in range(2):
                    v_num_lists[i].append(v_nums[i])
                t2 = time.time()
                elapsed_time = t2 - tb
                print(f"now {generation} generation, and spent {elapsed_time}")
                print(v_nums)
                tb = t2
            elif (v_num_lists_sorted[0][2] < v_nums[0] and
                  v_num_lists_sorted[0][7] > v_nums[0] and
                  v_num_lists_sorted[1][2] < v_nums[1] and
                  v_num_lists_sorted[1][7] > v_nums[1]):
                    print(f"simulation finish at {generation} generation")
                    focal = False
            else:
                t2 = time.time()
                elapsed_time = t2 - tb
                print(f"now {generation} generation, and spent {elapsed_time}")
                print(v_nums)
                for i in range(2):
                    del v_num_lists[i][0]
                    v_num_lists[i].append(v_nums[i])
                tb = t2
    t2 = time.time()
    elapsed_time = t2 - t1
    result = population.print_summary(parameter_obj)
    print(result)
    print(f"all time spent {elapsed_time}")
    return(*result, generation, elapsed_time)


file = ""
test_parameters = parameter_object(N=100000,
                                   mutation_rate_coef=1,
                                   mutater_effect=10,
                                   s_exp_mean=0.5,
                                   selection_coef=0.5,
                                   tsgnon_s=2,
                                   mutater_s=5)
with open(file, "w") as f:
    out = csv.writer(f, delimiter="\t")
    col = ["num", "TSG_nonsyn_num", "TSG_syn_num"]
    col += ["TSG_nonsyn_regression", "TSG_syn_regression"]
    col += ["generation", "spent_time"]
    out.writerow(col)
    for times in range(5):
        results = simulation(test_parameters)
        out.writerow([times, *results])
