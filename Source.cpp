#include <iostream>
#include <string>
#include <cmath>
#include <list>
#include <vector>
#include <valarray>

#include <chrono>

#include "NSim.h"

int main() {
    std::vector<class Particle> ps_tree, ps_direct;
    NSim_Init(ps_tree, 1000);
    NSim_Init(ps_direct, 1000);
    
    std::list<class Particle*> ptrs;
    std::list<class Particle*>::iterator left, right, part;

    double dt = 1 / 60.0;
    int i = 0;
    while (i<10) {
        auto start_clock = std::chrono::high_resolution_clock::now();
        
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
        
        auto end_clock = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_interval = end_clock - start_clock;
        std::cout << "T: " << time_interval.count() << "\t";

        start_clock = std::chrono::high_resolution_clock::now();
        NSim_Step(ps_direct, nullptr, dt);
        end_clock = std::chrono::high_resolution_clock::now();
        time_interval = end_clock - start_clock;
        std::cout << "D: " << time_interval.count() << std::endl;
	    //std::cout << ps_direct[0].pos[0] << std::endl;
        i++;
    }
}
