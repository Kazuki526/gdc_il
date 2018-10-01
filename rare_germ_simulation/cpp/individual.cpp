#include"individual.hpp"

Individual::Individual(int m, std::vector<int> tn){
  mutater=m;
  if(!tn.empty()){
    for(int &x:tn){
      if(count(tn.begin(), tn.end(), x) ==2){
        if(find(tsg_non_hom.begin(),tsg_non_hom.end(),x) != tsg_non_hom.end()){
          tsg_non_hom.push_back(x);
        }
      }else{
        tsg_non_het.push_back(x);
      }
    }
  }
}

void Individual::set_param(Parameters& param){
  mut_r = mutation_rate;
  mut_r *= param.mutation_rate_coef;
  for(int i=1; mutater >= i; i++){mut_r*=param.mutater_effect;}
  damage=0;
  for(int& mu:tsg_non_het){damage+=param.tsg_non_damage[mu-1];}
  for(int& mu:tsg_non_hom){damage+=param.tsg_non_damage[mu-1]*2;}
  fitness=1 - damage*param.fitness_coef;
}

void Individual::add_mutations(Parameters& param, std::mt19937& mt){
  // new mutater mutation
  std::bernoulli_distribution p_mutater(param.mutater_mutation_rate);
  int new_mutater=p_mutater(mt);
  if(new_mutater>2){new_mutater=2;}
  mutater+=new_mutater;
  if(mutater>2){
    mutater-=2;
    if(mutater==2){mutater=0;}
  }
  // new TSG nonsynonymous mutation
  std::poisson_distribution<> pois_tn(tsg_non_site*mut_r);
  std::uniform_int_distribution<> tn_mut(1, tsg_non_site);
  int tn_num=pois_tn(mt);
  std::vector<int> new_tsg_non ={};
  while(tn_num > 0){
    tn_num--;
    new_tsg_non.push_back(tn_mut(mt));
  }
  if(!new_tsg_non.empty()){
    for(int &mut:new_tsg_non){
      // mutation occur at homo site
      int focal=-1;
      for(int i=0; tsg_non_hom.size() > i; i++){
        if(tsg_non_hom[i] == mut){focal=i;}
      }
      if(focal != -1){
        tsg_non_hom.erase(tsg_non_hom.begin()+focal);
        tsg_non_het.push_back(mut);
        continue;
      }
      // mutation occur at hetero site
      // focal=-1;
      for(int mu=0; tsg_non_het.size() > mu; mu++){
        if(tsg_non_het[mu] == mut){focal=mu;}
      }
      if(focal != -1){
        tsg_non_het.erase(tsg_non_het.begin()+focal);
        tsg_non_hom.push_back(mut);
        continue;
      }
      tsg_non_het.push_back(mut);
    }
  }
}
