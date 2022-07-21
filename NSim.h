#ifndef NSIM_HH
#define NSIM_HH

//#include <valarray>
//#include <vector>
#include <random>
#include <iostream>

#include "Particle.h"
#include "MassTree.h"
#include "PoissonSolver.h"
#include "Integrator.h"

int Energy(std::vector<class Particle>& ps);
void NSim_Step(std::vector<class Particle>& ps, class Tree* T, double dt);
void NSim_Init(std::vector<class Particle>& ps, int num_ps);

#endif
