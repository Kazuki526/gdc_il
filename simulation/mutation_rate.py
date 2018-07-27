#! /usr/local/bin/env python
# -*- coding: utf-8 -*-
import random
import numpy as np

mutater_effect = 10
# defined parameters
TSG_nonsyn_site =
TSG_syn_site =
control_nonsyn_site =
control_syn_site =


def new_mutation(mp, site_num):
    new_mus = []
    for x in range(3):
        mutation_rate = 1.5*(10**-8) * mutater_effect ** x * site_num
        new_mus.extend(np.random.poisson(mutation_rate, mp[x]).tolist)
    return new_mus


class Individual:
    def __init__(self, mutater=0,
                 tsg_n=0, tsg_s=0, control_n=0, control_s=0):
        self._fitness = 1
        self._mutater = mutater
        self._tsg_n = tsg_n
        self._tsg_s = tsg_s
        self._control_n = control_n
        self._control_s = control_s

    @property
    def fitness(self):
        return self._fitness

    @property
    def mutater(self):
        return self._mutater

    # get [mutater, tsg_n, tsg_s, control_n, control_s]
    def mutations(self):
        mutations = [self._mutater, self._tsg_n, self._tsg_s]
        mutations += [self._control_n, self._control_s]
        return(mutations)

    # get haploid mutaion num list
    def reproduct(self):
        new_mutations = [np.random.binomial(mu, 0.5) for mu in self.mutations]
        return(new_mutations)


class Population:
    def __init__(self, population_N):
        for i in range(0, population_N):
            self.individuals.append(Individual())

    def get_fitness_list(self):
        fitness_list = []
        for i in len(self.individuals):
            fitness_list.append(i.fitness)
        return fitness_list

    def get_munum(self):
        muter_pnum = [len([x for x in self.individuals if x.mutater == 0])]
        muter_pnum.append(len([x for x in self.individuals if x.mutater == 1]))
        muter_pnum.append(len([x for x in self.individuals if x.mutater == 2]))
        new_mut_tn = new_mutation(muter_pnum, TSG_nonsyn_site)
        for i in range(len(new_mut_tn)):
            new_mus = random.choice(TSG_nonsyn_site, new_mut_tn[i])
            self.individuals[i]._tsg_n.append(new_mus)

    def next_generation_wf(self, population_N):
        fitnesses = self.get_fitness_list()

        next_generation = []
        for n in range(population_N):
            mutation1 = random.choice(self.individuals,
                                      weights=fitnesses).reproduct()
            mutation2 = random.choice(self.individuals,
                                      weights=fitnesses).reproduct()
            next_generation.append(Individual(*(mutation1 + mutation2)))
        self.individuals = next_generation
