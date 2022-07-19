#include "Particle.h"

Particle::Particle(double mass, std::valarray<double> coordinates) :
	m{ mass },
	coords{ coordinates },
	tmp_coords{ coordinates },
	pos{coordinates[pos_i]},
	vel{coordinates[vel_i]},
	pos_half{0,0,0},
	force{ 0,0,0 }
{}

void Particle::clearForce(void) {
	force = std::valarray<double>{ 0,0,0 };
}
