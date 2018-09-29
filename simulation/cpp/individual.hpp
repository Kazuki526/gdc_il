#pragma once
#ifndef IND_HPP
#define IND_HPP

#include<iostream>
#include<vector>
#include<random>
using namespace std;

extern const double mutation_rate;
extern const int tsg_non_site;
extern const int tsg_syn_site;
extern const int cont_non_site;

struct Parameters
{
  double mutation_rate_coef=10;
  double mutater_effect=100;
  double mutater_mutation_rate=0.001;
  double mutater_damage=0.3;
  double tsg_non_damage_e=0.5;
  double cont_non_damage_e=0.1;
  std::vector<double> tsg_non_damage;
  std::vector<double> cont_non_damage;
  bool cont_non_fitness_only=true;
  double fitness_coef=0.1;
  double cancer_prob_coef=1;
  void set_damage(){
    std::random_device rnd;
    std::mt19937 mt(rnd());
    std::exponential_distribution<> tsg_ex(tsg_non_damage_e);
    std::exponential_distribution<> cont_ex(cont_non_damage_e);
    for(int s=0; tsg_non_site > s; s++) {tsg_non_damage.push_back(tsg_ex(mt));}
    for(int s=0; cont_non_site > s; s++){cont_non_damage.push_back(cont_ex(mt));}
  }
};

class Individual
{
private:
  int mutater;
  std::vector<int> tsg_non_het;
  std::vector<int> tsg_non_hom;
  double mut_r;
  double damage;
  double fitness;
public:
  Individual():mutater(0),mut_r(mutation_rate),
               damage(0) ,fitness(1){};
  Individual(int m, std::vector<int> tn);
  int get_mutater(){return mutater;}
  std::vector<int> get_tsg_non_het(){return tsg_non_het;}
  std::vector<int> get_tsg_non_hom(){return tsg_non_hom;}
  double get_damage(){return damage;}
  double get_fitness(){return fitness;}
  void set_param(Parameters& param);
  void add_mutations(Parameters& param, std::mt19937& mt);

};

#endif
