#pragma once
#ifndef POP_HPP
#define POP_HPP

#include"individual.hpp"

extern const int N;
extern const double mutation_rate;
extern const int tsg_non_site;
extern const int tsg_syn_site;
extern const int cont_non_site;

class Population : public Individual
{
private:
  std::vector<Individual> individuals;
  std::mt19937 mt;
public:
  Population(Parameters& param);
  std::vector<int> get_mutater_list();
  void add_new_mutations(Parameters& param);
  Individual reproduct(Parameters& param,int i1, int i2);
  void next_generation(Parameters& param);
  std::vector<int> tsg_non_mutation_count();
};

#endif
