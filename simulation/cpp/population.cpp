#include<iostream>
#include<vector>
#include<random>
using namespace std;
const int N = 10;

class Individual
{
private:
  int id;
public:
  Individual(int i=1){id = i;}
  void set_id(int i){id = i;}
  int get_id(){return id;}
};

class Population : public Individual
{
private:
  Individual individuals[N];
public:
  Population();
  void set_id(int i, int id){individuals[i].set_id(id);}
  Individual get_ind(int i){return individuals[i];}
  void get_id_list(int ids[N]);
  void next_generation();
};

Population::Population(){
  for(int id=1; 10 >= id;id++){
    individuals[id].set_id(id);
  }
}

void Population::get_id_list(int ids[N]){
  for(int i=0; N > i; i++){
    ids[i] = individuals[i].get_id();
  }
}

void Population::next_generation(){
  Individual next_inds[N];
  std::random_device rnd;
  std::mt19937 mt(rnd());
  std::uniform_int_distribution<> dist(0, N -1);
  for(int i=0; N > i; i++){
     next_inds[i] = individuals[dist(mt)];
  }
  for(int i=0; N > i; i++){
    individuals[i] = next_inds[i];
  }
}


int main()
{
  Population population;
  population.next_generation();
  int id[N];
  population.get_id_list(id);
  for(int i=0; N > i; i++){
    cout <<id[i] <<" ";
  }
  cout <<endl;

  return 1;
}
