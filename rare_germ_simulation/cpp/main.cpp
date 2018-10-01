#include"individual.hpp"
#include"population.hpp"
using namespace std;

const int N = 50000;
const double mutation_rate = 0.00000001;
const int tsg_non_site = 191624;
const int tsg_syn_site = 61470;
const int cont_non_site = 923307;
const double mean_onset_age = 61.5;
const double age_sd = 13.5;

int main()
{
  Parameters param;
  param.set_damage();
  Population population(param=param);
  for(int t=1; 50 >= t; t++){
    population.add_new_mutations(param);
    population.next_generation(param);
    int rare_mut_num = 0;
    for(int &n:population.tsg_non_mutation_count()){
      if((n>0) && (n<(N*0.05/100))){rare_mut_num+=n;}
    }
    cout << rare_mut_num << endl;
  }
  return 0;
}
