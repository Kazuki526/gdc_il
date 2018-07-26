#! /usr/local/bin/env python
# -*- coding: utf-8 -*-
import random
import numpy as np


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
        self.individuals = Individual(mutater=1)
        for i in range(1, population_N):
            self.individuals.append(Individual())

    def get_fitness_list(self):
        fitness_list = []
        for i in len(self.individuals):
            fitness_list.append(i.fitness)
        return fitness_list

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
