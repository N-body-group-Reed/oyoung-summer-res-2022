#include "PoissonSolver.h"

double softeningLength;

template<typename T>
T innerProduct(const std::valarray<T>& v1, const std::valarray<T>& v2) {
	return (v1 * v2).sum();
}

void Poisson_DirectSummation(std::vector<class Particle>& ps) {
	for (std::vector<class Particle>::iterator i = ps.begin(); i != ps.end(); i++) {
		for (std::vector<class Particle>::iterator j = i + 1; j != ps.end(); j++) {
			Poisson_CalculateForce_Newton(*i,*j);
		}
	}
}

void Poisson_CalculateForce_Hook(class Particle& p1, class Particle& p2) {
	std::valarray<double> r = std::valarray<double>(p1.tmp_coords[pos_i]) - std::valarray<double>(p2.tmp_coords[pos_i]);
	std::valarray<double> F = -r;
	p1.force += F;
	p2.force += -F;
}

void Poisson_CalculateForce_Newton(class Particle& p1, class Particle& p2) {
	std::valarray<double> r = std::valarray<double>(p1.tmp_coords[pos_i]) - std::valarray<double>(p2.tmp_coords[pos_i]);
	double R = pow(sqrt(innerProduct(r, r)), 3);
	std::valarray<double> F = -G * p1.m * p2.m * r / R;
	p1.force += F;
	p2.force += -F;
}

void Poisson_CalculateForce_SoftSphere(class Particle& p1, class Particle& p2) {
	std::valarray<double> r = std::valarray<double>(p1.tmp_coords[pos_i]) - std::valarray<double>(p2.tmp_coords[pos_i]);
	double SF = (innerProduct(r, r) + 5 * pow(softeningLength, 2) / 2) / pow(sqrt(innerProduct(r, r) + pow(softeningLength, 2)), 5);
	std::valarray<double> F = -G * p1.m * p2.m * r * SF;
	p1.force += F;
	p2.force += -F;
}

void Poisson_TreeForce(std::vector<class Particle>& ps, class Tree* T) {
	for (std::vector<class Particle>::iterator p = ps.begin(); p != ps.end(); p++) {
		Poisson_CalculateForce_Tree(*p, T);
	}
}

const double theta = 0.75;
void Poisson_CalculateForce_Tree(class Particle& p, class Tree* T) {
	std::valarray<double> r = p.pos - T->r_cm;
	//double SF = (innerProduct(r, r) + 5 * pow(softeningLength, 2) / 2) / pow(sqrt(innerProduct(r, r) + pow(softeningLength, 2)), 5);
	double com_dist = sqrt(innerProduct(r, r));
	//std::cout << com_dist << std::endl;
	if (T->l != nullptr && T->r != nullptr) {
		//std::cout << "Has Children" << std::endl;
		// cell length would be better quantified by the diagonal length of box
		// more difficult to calculate with max difference approx 1.4
		// how good of a measure is "node length" anyway?
		double cell_length;
		if (T->axis == Partition_axis::X) cell_length = T->x_max - T->x_min;
		else if (T->axis == Partition_axis::Y) cell_length = T->y_max - T->y_min;
		else cell_length = T->z_max - T->z_min;

		if (cell_length / com_dist < theta) {
			//std::cout << "N*ln N interaction" << std::endl;
			double R = pow(com_dist, 3);
			p.force += -G * p.m * T->mass * r / R; //* SF;
		}
		else {
			//std::cout << "open children" << std::endl;
			Poisson_CalculateForce_Tree(p, T->l);
			Poisson_CalculateForce_Tree(p, T->r);
		}
	}
	else {
		double tol = 0.1;
		if (com_dist > tol) {
			//std::cout << "N*N interaction" << std::endl;
			double R = pow(com_dist, 3);
			p.force += -G * p.m * T->mass * r / R; //* SF;
		}
	}
}

void PoissonSolver(std::vector<class Particle>& ps, class Tree* T) {
	softeningLength = .98 * pow(ps.size(), -0.26);
	if (T == nullptr) Poisson_DirectSummation(ps);
	else Poisson_TreeForce(ps, T);
}