#ifndef PARTICLE_HH
#define PARTICLE_HH

#include <valarray>

// RK4
const std::slice pos_i = std::slice(0, 3, 1);
const std::slice vel_i = std::slice(3, 3, 1);

class Particle {
public:
	int tag;
	double m;
	// RK4
	std::valarray<double> coords;
	std::valarray<double> tmp_coords;
	std::valarray<double> k1;
	std::valarray<double> k2;
	std::valarray<double> k3;
	std::valarray<double> k4;
	// LeapFrog
	std::valarray<double> pos;
	std::valarray<double> vel;
	std::valarray<double> pos_half;
	std::valarray<double> force;

	Particle(double /*mass*/ = 1.0, std::valarray<double> /*position,velocity*/ = std::valarray<double>{ 0,0,0, 0,0,0 });
	void clearForce(void);
};

#endif