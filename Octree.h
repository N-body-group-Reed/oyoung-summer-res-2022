#ifndef OCTREE_HH
#define OCTREE_HH

#include <list>
#include <vector>
#include <iostream>

#include "Particle.h"

class OctTree {
	class TreeNode* root;
};

class TreeNode {
public:
	double mass;
	std::valarray<double> r_cm;
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;
	
	class Tree* l;
	class Tree* r;
	Tree(std::list<class Particle*>& ps, std::list<class Particle*>::iterator& start, std::list<class Particle*>::iterator& stop);
	~Tree(void);
	void computeBondingBox(std::list<class Particle*>::iterator start, std::list<class Particle*>::iterator stop);
	void partitionVolume(void);
	void sortParticles(std::list<class Particle*>& ps, std::list<class Particle*>::iterator& start, std::list<class Particle*>::iterator& stop);
	std::list<class Particle*>::iterator getPartitionIterator(std::list<class Particle*>::iterator start, std::list<class Particle*>::iterator stop);
	void computeMassMoments(std::list<class Particle*>::iterator start, std::list<class Particle*>::iterator stop);
	// center of mass and mass moments
};

void buildTree(class Tree* T, std::list<class Particle*>& ptrs, std::list<class Particle*>::iterator start, std::list<class Particle*>::iterator stop, std::list<class Particle*>::iterator part);
void orderParticles(std::vector<class Particle>& ps, std::list<class Particle*>& ptrs);
void printLeafs(class Tree* T);

#endif
