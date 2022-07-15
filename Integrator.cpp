#include "Integrator.h"

void Integrator_RK4(std::vector<class Particle>& ps, double dt) {
	PoissonSolver(ps);
	for (std::vector<class Particle>::iterator p = ps.begin(); p != ps.end(); p++) {
		p->k1 = dt * RK_Vector(*p);
		p->clearForce();
		p->tmp_coords += p->k1 / 2;
	}
	PoissonSolver(ps);
	for (std::vector<class Particle>::iterator p = ps.begin(); p != ps.end(); p++) {
		p->k2 = dt * RK_Vector(*p);
		p->clearForce();
		p->tmp_coords = p->coords + p->k2 / 2;
	}
	PoissonSolver(ps);
	for (std::vector<class Particle>::iterator p = ps.begin(); p != ps.end(); p++) {
		p->k3 = dt * RK_Vector(*p);
		p->clearForce();
		p->tmp_coords = p->coords + p->k3;
	}
	PoissonSolver(ps);
	for (std::vector<class Particle>::iterator p = ps.begin(); p != ps.end(); p++) {
		p->k4 = dt * RK_Vector(*p);
		p->clearForce();
		p->tmp_coords = p->coords + (1.0 / 3.0) * (p->k1 / 2 + p->k2 + p->k3 + p->k4 / 2);
		p->coords = p->tmp_coords;
	}
}

std::valarray<double> RK_Vector(class Particle& p) {
	std::valarray<double> vec(6);
	vec[pos_i] = p.tmp_coords[vel_i];
	vec[vel_i] = p.force / p.m;
	return vec;
}

void Integrator_LeapFrog_DKD(std::vector<class Particle>& ps, class Tree* T, double dt) {
	for (std::vector<class Particle>::iterator p = ps.begin(); p != ps.end(); p++) {
		p->pos_half = p->pos + dt * p->vel / 2;
		// artifact of PoissonSolver for RK4
		std::valarray<double> tmp(6);
		tmp[pos_i] = p->pos_half;
		tmp[vel_i] = p->vel;
		p->tmp_coords = tmp;
	}
	PoissonSolver(ps,T);
	for (std::vector<class Particle>::iterator p = ps.begin(); p != ps.end(); p++) {
		p->vel += dt * p->force / p->m;
		p->pos = p->pos_half + dt * p->vel / 2;
		p->clearForce();
		// artifact of PoissonSolver for RK4
		std::valarray<double> tmp(6);
		tmp[pos_i] = p->pos;
		tmp[vel_i] = p->vel;
		p->coords = tmp;
	}
}


void Integrator(std::vector<class Particle>& ps, class Tree* T, double dt) {
	Integrator_LeapFrog_DKD(ps,T,dt);
}