#include <iostream>
#include <string>
#include <cmath>
#include <list>
#include <vector>
#include <valarray>

#include "NSim.h"

int main() {
    std::vector<class Particle> ps_tree, ps_direct;
    NSim_Init(ps_tree, 15);
    NSim_Init(ps_direct, 15);
    
    std::list<class Particle*> ptrs;
    std::list<class Particle*>::iterator left, right, part;

    double dt = 1 / 60.0;
    int i = 0;
    while (i<1000) {
        for (std::vector<class Particle>::iterator p = ps_tree.begin(); p != ps_tree.end(); p++) ptrs.push_back(&*p);
        left = ptrs.begin();
        right = ptrs.end();
        Tree* T = new Tree(ptrs, left, right);
        part = T->getPartitionIterator(left, right);
        buildTree(T, ptrs, left, right, part);
        T->computeMassMoments(left, right);

        orderParticles(ps_tree, ptrs);

        NSim_Step(ps_tree, T, dt);

        delete T;
        ptrs.erase(ptrs.begin(), ptrs.end());

        NSim_Step(ps_direct, nullptr, dt);
	std::cout << ps_direct[0].pos[0] << std::endl;
        i++;
    }
}
