#include"individual.hpp"
#include"population.hpp"

Population::Population(Parameters& param){
  for(int i=0; N > i; i++){
    Individual ind;
    ind.set_param(param);
    individuals.push_back(ind);
  }
  std::random_device rnd;
  std::mt19937 mt_(rnd());
  mt = mt_;
}

void Population::add_new_mutations(Parameters& param){
  for(int i=0; N > i; i++){
    individuals[i].add_mutations(param=param, mt=mt);
  }
}

Individual Population::reproduct(Parameters& param,int i1, int i2){
  int mutater=0;
  std::vector<int> tsg_non;
  std::bernoulli_distribution bern(0.5);
  // individuals[i1] 's gamate
  int muta = individuals[i1].get_mutater();
  if(muta == 2){
    mutater++;
  }else if(muta == 1){
    if(bern(mt)){mutater++;}
  }
  for(int &het_mu:individuals[i1].get_tsg_non_het()){
    if(bern(mt)){tsg_non.push_back(het_mu);}
  }
  for(int &hom_mu:individuals[i1].get_tsg_non_hom()){
    tsg_non.push_back(hom_mu);
  }
  // individuals[i2] 's gamate
  muta = individuals[i2].get_mutater();
  if(muta == 2){
    mutater++;
  }else if(muta == 1){
    if(bern(mt)){mutater++;}
  }
  for(int &het_mu:individuals[i2].get_tsg_non_het()){
    if(bern(mt)){tsg_non.push_back(het_mu);}
  }
  for(int &hom_mu:individuals[i2].get_tsg_non_hom()){
    tsg_non.push_back(hom_mu);
  }
  Individual ind(mutater, tsg_non); //明示的代入がダメ？
  ind.set_param(param);
  return(ind);
}

void Population::next_generation(Parameters& param){
  std::vector<Individual> next_inds;
  std::uniform_int_distribution<> dist(0, N -1);
  for(int i=0; N > i; i++){
     next_inds.push_back(reproduct(param=param, dist(mt), dist(mt))); //ここも明示的にできない
  }
  individuals = next_inds;
}

std::vector<int> Population::tsg_non_mutation_count(){
  std::vector<int> tsg_non_mutation(tsg_non_site);
  for(Individual &ind:individuals){
    for(int &het_mu:ind.get_tsg_non_het()){
      tsg_non_mutation[het_mu]++;
    }
    for(int &hom_mu:ind.get_tsg_non_hom()){
      tsg_non_mutation[hom_mu]+=2;
    }
  }
  return(tsg_non_mutation);
}
