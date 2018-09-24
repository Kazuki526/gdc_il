#include<iostream>
#include<vector>
using namespace std;
count int N = 10;

class Individual
{
private:
  int id;
public:
  Individual(int i=1){id = i;}
  void set_id(int i){id = i;}
}

class Population : public Individual
{
private:
  Individual individual[N]
public:
  Population()
  void next_generation()
}

Population::Population(){
  for(int id=1; 10 >= id;id++){
    individual[i] = Individual.set_id(i);
  }
}

void Population::next_generation(){
  
}

int main()
{
  Population population()

}
