#! /usr/local/bin/env python
# -*- coding: utf-8 -*-
import random
import numpy as np
import time
import copy
from sklearn import linear_model
from statistics import mean, stdev
from math import log10

# defined parameters
tsg_non_site = 191624
tsg_syn_site = 61470
cont_non_site = 923307
cont_syn_site = 319944
mean_onset_age = 61.5
age_sd = 13.5


# define class of parameters
class parameter_object:
    def __init__(self,
                 N,
                 mutation_rate_coef,
                 mutater_effect,
                 mutater_mutation_rate,
                 mutater_damage,
                 tsg_non_damage,
                 cont_non_damage,
                 cont_non_fitness_only,
                 fitness_coef,
                 cancer_prob_coef,
                 syn_damage=0):
        self.N = N
        self.mutation_rate = 1*(10**-8)*mutation_rate_coef
        self.mutater_effect = mutater_effect
        self.mutater_mutation_rate = mutater_mutation_rate
        tsgnon_d_list = np.random.exponential(tsg_non_damage, tsg_non_site)
        self.tsgnon_d = tsgnon_d_list.tolist()
        contnon_d_list = np.random.exponential(cont_non_damage, cont_non_site)
        self.contnon_d = contnon_d_list.tolist()
        syn_d_list = np.random.exponential(syn_damage, tsg_syn_site)
        self.syn_d = syn_d_list.tolist()
        self.cont_non_fitness_only = cont_non_fitness_only
        self.fitness_coef = fitness_coef
        self.cancer_prob_coef = cancer_prob_coef


def genotype_divide(mutations):
    homo = [x for x in set(mutations) if mutations.count(x) > 1]
    hetero = list(set(mutations) - set(homo))
    return(homo, hetero)


class Individual:
    def __init__(self, mutater, tsg_non, tsg_syn, cont_non, parameters):
        self._mutater = mutater
        self._tsg_non_hom, self._tsg_non_het = genotype_divide(tsg_non)
        self._tsg_syn_hom, self._tsg_syn_het = genotype_divide(tsg_syn)
        self._cont_non_hom, self._cont_non_het = genotype_divide(cont_non)
        damage = sum([parameters.tsgnon_d[mut] for mut in tsg_non])
        damage += sum([parameters.syn_d[mut] for mut in tsg_syn])
        if parameters.cont_non_fitness_only == "T":
            damage += sum([parameters.contnon_d[mut] for mut in cont_non])
        self._damage = damage
        fitness = 1 - damage * parameters.fitness_coef
        if parameters.cont_non_fitness_only == "F":
            cn_dam = sum([parameters.contnon_d[mut] for mut in cont_non])
            fitness -= cn_dam * parameters.fitness_coef
        self._fitness = 0 if fitness < 0 else fitness

    @property
    def mutater(self):
        return self._mutater

    @property
    def damage(self):
        return self._damage

    @property
    def fitness(self):
        return self._fitness

    # get [tsg_non_het, tsg_syn_het, cont_non_het, cont_syn_het]
    def mutations(self):
        mutations = copy.deepcopy([self._tsg_non_het, self._tsg_syn_het])
        mutations += copy.deepcopy([self._cont_non_het])
        return(mutations)

    # get [tsg_non_homo, tsg_syn_homo, cont_non_homo, cont_syn_homo]
    def homo_mutations(self):
        mutations = copy.deepcopy([self._tsg_non_hom, self._tsg_syn_hom])
        mutations += copy.deepcopy([self._cont_non_hom])
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

    def add_cont_non(self, muts):
        old_het = set(copy.copy(self._cont_non_het))
        old_hom = set(copy.copy(self._cont_non_hom))
        self._cont_non_het = list(old_het ^ set(muts))
        self._cont_non_het += list(old_hom & set(muts))
        self._cont_non_hom = list(old_hom ^ set(muts))
        self._cont_non_hom += list(old_het & set(muts))

    def add_mutater(self, add_num):
        bef = copy.copy(self._mutater)
        if bef+add_num != 2:
            self._mutater = (bef+add_num) % 2
        elif bef == 1 and add_num == 1:
            self._mutater = np.random.randint(2)
        else:
            self.mutater = 2

    def variant_num_tsg_non(self, variant_list):
        num = len(set(self._tsg_non_het) & set(variant_list))
        num += len(set(self._tsg_non_hom) & set(variant_list)) * 2
        return(num)

    def variant_num_tsg_syn(self, variant_list):
        num = len(set(self._tsg_syn_het) & set(variant_list))
        num += len(set(self._tsg_syn_hom) & set(variant_list)) * 2
        return(num)

    def variant_num_cont_non(self, variant_list):
        num = len(set(self._cont_non_het) & set(variant_list))
        num += len(set(self._cont_non_hom) & set(variant_list)) * 2
        return(num)


# make de nove mutations list (list of each ind de novo mutation num)
def new_mutation(mutater_num, site_num, params):
    mutation_rate = params.mutation_rate * site_num
    mut_ef = params.mutater_effect
    mut0 = np.random.poisson(mutation_rate, mutater_num.count(0)).tolist()
    mut1 = np.random.poisson(mutation_rate * mut_ef,
                             mutater_num.count(1)).tolist()
    mut2 = np.random.poisson(mutation_rate * (mut_ef**2),
                             mutater_num.count(2)).tolist()
    new_mus = []
    for i in range(params.N):
        if mutater_num[i] == 0:
            new_mus.append(mut0.pop())
        elif mutater_num[i] == 1:
            new_mus.append(mut1.pop())
        else:
            new_mus.append(mut2.pop())
    return new_mus


# make offspring from two Individuals
def reproduct(ind1, ind2, params):
    mutater = np.random.binomial(ind1.mutater + ind2.mutater, 0.5)
    mutater = 2 if mutater > 2 else mutater
    muts = [ind1.mutations()[i] + ind2.mutations()[i] for i in range(3)]
    new_mut = [random.sample(l, np.random.binomial(len(l), 0.5)) for l in muts]
    new_mut = [new_mut[i] + ind1.homo_mutations()[i] for i in range(3)]
    return(Individual(mutater, *new_mut, params))


class Population:
    def __init__(self, params):
        self.individuals = [Individual(mutater=0, tsg_non=[], tsg_syn=[],
                                       cont_non=[], parameters=params)
                            for i in range(params.N)]

    def get_fitness_list(self, params):
        fitness_list = [ind.fitness for ind in self.individuals]
        fit_sum = sum(fitness_list)
        fitness_list = [fit / fit_sum for fit in fitness_list]
        return fitness_list

    def get_cancer_prob(self, params):
        damage_list = [ind.damage for ind in self.individuals]
        prob_list = [1 + d * params.cancer_prob_coef for d in damage_list]
        prob_list = [0 if p < 0 else p for p in prob_list]
        prob_list = [prob / sum(prob_list) for prob in prob_list]
        return prob_list

    def get_mutater(self):
        mutater_list = [ind.mutater for ind in self.individuals]
        return(mutater_list)

    # add new mutations to each individuals
    def add_new_mutation(self, params):
        # individuals num of [mutater=0, mutater=1, mutater=2]
        mutater_num = self.get_mutater()
        new_mutater = np.random.binomial(1, params.mutater_mutation_rate,
                                         params.N).tolist()
        new_mut_tn = new_mutation(mutater_num, tsg_non_site, params)
        new_mut_ts = new_mutation(mutater_num, tsg_syn_site, params)
        new_mut_cn = new_mutation(mutater_num, cont_non_site, params)
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
            if new_mut_cn[n] != 0:
                new_mus = np.random.randint(0, cont_non_site,
                                            new_mut_cn[n]).tolist()
                self.individuals[n].add_cont_non(new_mus)

    # make next generation population
    def next_generation_wf(self, params):
        fitness = self.get_fitness_list(params=params)
        rand_ind = np.random.choice(self.individuals, params.N*2, p=fitness)
        next_generation = [reproduct(rand_ind[n], rand_ind[n+1], params)
                           for n in range(0, params.N*2, 2)]
        self.individuals = next_generation

    def variant_allelecount(self):
        v_AC = [[0 for i in range(tsg_non_site)]]
        v_AC += [[0 for i in range(tsg_syn_site)]]
        v_AC += [[0 for i in range(cont_non_site)]]
        for ind in self.individuals:
            het_muts = ind.mutations()
            hom_muts = ind.homo_mutations()
            for i in range(3):
                for het_mut in het_muts[i]:
                    v_AC[i][het_mut] += 1
                for hom_mut in hom_muts[i]:
                    v_AC[i][hom_mut] += 2
        return(v_AC)

    def print_rare_variant_num(self, params):
        v_AC = self.variant_allelecount()
        rare_nums = [0, 0, 0]
        for i in range(3):
            rare_nums[i] = sum([v_num for v_num in v_AC[i]
                                if v_num < params.N*2*0.0005])
        return(rare_nums)

    def print_summary(self, params, sample_num=6000):
        sample_num = params.N if sample_num > params.N else sample_num
        v_AC = self.variant_allelecount()
        rare_variants = [[], [], []]
        rare_nums = []
        for i in range(3):
            rare_variants[i] = [v for v in range(len(v_AC[i]))
                                if v_AC[i][v] < params.N*2*0.0005]
            rare_nums += [sum([v_AC[i][v] for v in rare_variants[i]])]
        cancer_prob = self.get_cancer_prob(params)
        sample = np.random.choice(self.individuals, sample_num,
                                  p=cancer_prob, replace=False).tolist()
        sample_age = [mean_onset_age-ind.damage for ind in sample]
        sample_age = [0 if age < 0 else age for age in sample_age]
        sample_tn_num = [ind.variant_num_tsg_non(rare_variants[0])
                         for ind in sample]
        sample_ts_num = [ind.variant_num_tsg_syn(rare_variants[1])
                         for ind in sample]
        sample_cn_num = [ind.variant_num_cont_non(rare_variants[2])
                         for ind in sample]
        sample_age = np.array(sample_age).reshape(-1, 1)
        sample_tn_num = np.array(sample_tn_num).reshape(-1, 1)
        sample_ts_num = np.array(sample_ts_num).reshape(-1, 1)
        sample_cn_num = np.array(sample_cn_num).reshape(-1, 1)
        tsg_non_reg = linear_model.LinearRegression()
        tsg_non_reg.fit(sample_age, sample_tn_num)
        tsg_syn_reg = linear_model.LinearRegression()
        tsg_syn_reg.fit(sample_age, sample_ts_num)
        cont_non_reg = linear_model.LinearRegression()
        cont_non_reg.fit(sample_age, sample_cn_num)
        return(*rare_nums,
               float(tsg_non_reg.coef_), float(tsg_syn_reg.coef_),
               float(cont_non_reg.coef_))


def simulation(parameter_obj):
    t1 = time.time()
    population = Population(params=parameter_obj)
    focal = True
    tsgnon_v_num = []
    tsgsyn_v_num = []
    contnon_v_num = []
    generation = 0
    while focal:
        generation += 1
        population.add_new_mutation(parameter_obj)
        population.next_generation_wf(parameter_obj)
        v_nums = population.print_rare_variant_num(parameter_obj)
        if generation <= 100:
            tsgnon_v_num.append(v_nums[0])
            tsgsyn_v_num.append(v_nums[1])
            contnon_v_num.append(v_nums[2])
        else:
            del tsgnon_v_num[0]
            del tsgsyn_v_num[0]
            del contnon_v_num[0]
            tsgnon_v_num.append(v_nums[0])
            tsgsyn_v_num.append(v_nums[1])
            contnon_v_num.append(v_nums[2])
            tn = stdev(tsgnon_v_num)/mean(tsgnon_v_num)
            ts = stdev(tsgsyn_v_num)/mean(tsgsyn_v_num)
            cn = stdev(contnon_v_num)/mean(contnon_v_num)
            if (tn < 0.05 and ts < 0.05 and cn < 0.05) or generation > 500:
                focal = False
        # if generation % 10 == 0:
        #     t2 = time.time()
        #     elapsed_time = t2 - t1
        #     muter = population.get_mutater()
        #     muter_num = [muter.count(0), muter.count(1), muter.count(2)]
        #     print(f"now {generation} : {elapsed_time} sec:")
        #     print(v_nums, muter_num)
    result = population.print_summary(parameter_obj)
    t2 = time.time()
    elapsed_time = t2 - t1
    print(f"all time spent {generation} : {elapsed_time} sec")
    return(result)


# test_parameters = parameter_object(N=10000,
#                                    mutation_rate_coef=10,
#                                    mutater_effect=50,
#                                    mutater_mutation_rate=1*(10**-3),
#                                    mutater_damage=5,
#                                    tsg_non_damage=1,
#                                    cont_non_damage=1,
#                                    cont_non_fitness_only="F",
#                                    fitness_coef=0.1,
#                                    cancer_prob_coef=0.001,
#                                    syn_damage=0)
def log_rand(min, max):
    log_min = log10(min) if min != 0 else log10(0.000001)
    log_max = log10(max)
    x = random.uniform(log_min, log_max)
    x = 10 ** x if x != 0.000001 else 0
    return(x)


def rand_params():
    rand = random.uniform
    p_list = [50000,                      # N
              rand(0.1, 1),                # mutation_rate_coef
              rand(5, 100),                # mutation_effect
              rand(0, 0.001),              # mutater_mutation_rate
              rand(0, 5),                 # mutater_damage
              rand(0.1, 2),               # tsg_non_damage
              rand(0.1, 1),               # cont_non_damage
              random.choice(["T", "F"]),  # cont_non_fitness_only
              log_rand(0.001, 0.5),       # fitness_coef
              log_rand(0, 0.1)]           # cancer_prob_coef
    return(p_list)


with open("pre_run.tsv", "w", 1) as f:
    cols = ["N", "mutation_rate_coef", "mutater_effect", "mutater_damage",
            "tsg_non_damage", "cont_non_damage", "cont_non_fitness_only",
            "fitness_coef", "cancer_prob_coef",
            "tsg_non_num", "tsg_syn_num", "cont_non_num",
            "tsg_non_reg", "tsg_syn_reg", "cont_non_reg"]
    print(*cols, sep='\t', file=f)
    for times in range(100):
        print(f"simulate {times}")
        rand_parameters = rand_params()
        out = simulation(parameter_object(*rand_parameters))
        print(*rand_parameters, *out, sep='\t', file=f)
