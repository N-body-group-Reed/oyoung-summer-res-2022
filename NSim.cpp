#include "NSim.h"

double EnergyInit;

int Energy(std::vector<class Particle>& ps) {
  double TE = 0;
  for (std::vector<class Particle>::iterator p = ps.begin(); p != ps.end(); p++) {
    std::valarray<double> v = p->vel;
    double KE =1.0/2.0*innerProduct(v, v)*p->m;
    //std::cout<<"KE"<<KE<<std::endl;
    double UE=0;
    for (std::vector<class Particle>::iterator p1 = ps.begin(); p1 != ps.end(); p1++) {
      if (p != p1){
	std::valarray<double> r = p->pos - p1->pos;
	double R = sqrt(innerProduct(r, r));
	UE += -G * p->m * p1->m /R ;
      }
    }
    //std::cout<<"UE"<<UE<<std::endl;
    //std::cout<<"E"<<UE+KE<<std::endl;
    TE+=UE+KE;
  }
  return TE;
}
std::valarray<double> momentum(std::vector<class Particle>& ps) {
  std::valarray<double> TP = {0.0,0.0,0.0};
  for (std::vector<class Particle>::iterator p = ps.begin(); p != ps.end(); p++) {
    std::valarray<double> v = p->vel;
    std::valarray<double> P =v * p->m;
  TP+=P;
  }
  return TP;
}

std::valarray<double> spin(std::vector<class Particle>& ps) {
  std::valarray<double> TP = {0.0,0.0,0.0};
  std::valarray<double> com = {0.0,0.0,0.0};
  double totalmass=0;
  for (std::vector<class Particle>::iterator p = ps.begin(); p != ps.end(); p++) {
    totalmass+=p->m;
    com+=p->pos*p->m;
  }
  com=com/totalmass;
  for (std::vector<class Particle>::iterator p = ps.begin(); p != ps.end(); p++) {
    std::valarray<double> v = p->vel;
    std::valarray<double> r =p->pos-com;
    std::valarray<double> P ={v[2]*r[3]-v[3]*r[2],   v[3]*r[1]-v[1]*r[3],   v[1]*r[2]-v[2]*r[1]};
    P=p->m*P;
  TP+=P;
  }
  return TP;
}
double mag(std::valarray<double> P){
  double mag;
  for (int i=0; i != 3;i++) mag+=pow(P[i],2);
  mag=pow(mag,0.5);
  return mag;
}

void NSim_Step(std::vector<class Particle>& ps, class OCTree* T, double dt) {
	Integrator(ps, T, dt);
	if (T != nullptr) std::cout<<(Energy(ps)-EnergyInit)/EnergyInit<<std::endl;
	std::cout<<"momentum "<<mag(momentum(ps)-momentumInit)/mag(momentumInit)<<std::endl;
	std::cout<<"spin "<<mag(spin(ps)-spinInit)/mag(spinInit)<<std::endl;
}

void NSim_Init(std::vector<class Particle>& ps, int num_ps) {
    std::default_random_engine generator;
    std::exponential_distribution<double> mass_dist(.3);
    std::uniform_real_distribution<double> pos_dist(-10.0, 10.0);
    std::uniform_real_distribution<double> vel_dist(-2.5, 2.5);
    std::vector<class Particle> tmp;
    for (int i = 0; i < num_ps; i++) {
        double mass = mass_dist(generator) + 1;
        //std::cout << mass << std::endl;
        std::valarray<double> coords{ pos_dist(generator),pos_dist(generator),pos_dist(generator),vel_dist(generator),vel_dist(generator),vel_dist(generator) };
        Particle p(mass, coords);
        p.tag = i;
        tmp.push_back(p);
    }
	ps.assign(tmp.begin(), tmp.end());
	EnergyInit = Energy(ps);	
}
