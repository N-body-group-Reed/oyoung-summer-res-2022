#ifndef INTEGRATOR_HH
#define INTEGRATOR_HH

#include "Particle.h"
#include "PoissonSolver.h"

void Integrator_RK4(std::vector<class Particle>& ps, double dt);
std::valarray<double> RK_Vector(class Particle& p);

void Integrator_LeapFrog_DKD(std::vector<class Particle>& ps, class Tree* T, double dt);

//void Integrator_LeapFrog_KDK(std::vector<class Particle>& ps, double dt);

void Integrator(std::vector<class Particle>& ps, class Tree* T, double dt);

#endif