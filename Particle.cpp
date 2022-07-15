#include "Particle.h"

Particle::Particle(double mass, std::valarray<double> coordinates) :
	m{ mass },
	coords{ coordinates },
	force{ 0,0,0 },
	tmp_coords{ coordinates },
	pos_half{0,0,0},
	pos{coordinates[pos_i]},
	vel{coordinates[vel_i]}
{}

void Particle::clearForce(void) {
	force = std::valarray<double>{ 0,0,0 };
}