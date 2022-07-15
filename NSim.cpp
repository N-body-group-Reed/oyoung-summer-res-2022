#include "NSim.h"

void NSim_Step(std::vector<class Particle>& ps, class Tree* T, double dt) {
	Integrator(ps, T, dt);
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
}