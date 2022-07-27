#include "MassTree.h"

Tree::Tree(std::list<class Particle*>& ps, std::list<class Particle*>::iterator& start, std::list<class Particle*>::iterator& stop,int index,int depth) :
	mass{0},
	r_cm{0,0,0},
	x_min{ (*start)->pos[0] },
	x_max{ (*start)->pos[0] },
	y_min{ (*start)->pos[1] },
	y_max{ (*start)->pos[1] },
	z_min{ (*start)->pos[2] },
	z_max{ (*start)->pos[2] },
	l{nullptr},
	r{nullptr}
{
}

Tree::~Tree(void) {
	if (l != nullptr) delete l;
	if (r != nullptr) delete r;
}

void Tree::computeBondingBox(std::list<class Particle*>::iterator start, std::list<class Particle*>::iterator stop) {
	//std::cout << "new tree" << std::endl;
	for (std::list<class Particle*>::iterator p = start; p != stop; p++) {
		if ((*p)->pos[0] < x_min) x_min = (*p)->pos[0];
		if ((*p)->pos[0] > x_max) x_max = (*p)->pos[0];
		if ((*p)->pos[1] < y_min) y_min = (*p)->pos[1];
		if ((*p)->pos[1] > y_max) y_max = (*p)->pos[1];
		if ((*p)->pos[2] < z_min) z_min = (*p)->pos[2];
		if ((*p)->pos[2] > z_max) z_max = (*p)->pos[2];
	}
	//std::cout << "x: (" << x_min << "," << x_max << ")\t";
	//std::cout << "y: (" << y_min << "," << y_max << ")\t";
	//std::cout << "z: (" << z_min << "," << z_max << ")" << std::endl;
}



int Tree::BinaryDigit(double N,int D){
  double n=N/pow(2,-D+1);
  return (int)(2*(n-(int)(n)));
}

double Tree::Normalize(double x,int index){
  double zero = pow(10,-5);
  if (index==0) x=(x-x_min)/(x_max-x_min)+zero*x;
  else if (index==1) x=(x-y_min)/(y_max-y_min)+zero*x;
  else x=(x-z_min)/(z_max-z_min)+zero*x;
  return x;
}

std::list<class Particle*>::iterator Tree::sortParticles(std::list<class Particle*>& ps, std::list<class Particle*>::iterator& start, std::list<class Particle*>::iterator& stop,int index,int depth) {
  computeBondingBox(start,stop);
  int dist;
  int zeroes=0;
  for (std::list<class Particle*>::iterator itr = start; itr != stop; itr++) {
    //std::cout<<Normalize((**itr).pos[index],index)<<",";
    //std::cout<<(BinaryDigit(Normalize((**itr).pos[index],index),depth+1)==0)<<",";
    if (BinaryDigit(Normalize((**itr).pos[index],index),depth+1)==0){
      dist=std::distance(start,itr);
      ps.splice(start, ps, itr);
      start=itr;
      std::advance(itr,dist);
      zeroes+=1;
    }
  }
  std::list<class Particle*>::iterator part = start;
  std::advance(part,zeroes-1);
  //std::cout<<"sorted"<<distance(start,stop)<<",";
  //std::cout<<"sorted"<<distance(part,stop)<<",";
  //std::cout<<"sorted"<<distance(start,part)<<std::endl;
  return part;
}



void Tree::computeMassMoments(std::list<class Particle*>::iterator start, std::list<class Particle*>::iterator stop) {
	if (l == nullptr && r == nullptr) {
		for (std::list<class Particle*>::iterator p = start; p != stop; p++) {
			mass += (*p)->m;
			r_cm += (*p)->m * (*p)->pos;
		}
		r_cm = r_cm / mass;
	}
	else if (l != nullptr && r != nullptr) {
		mass = l->mass + r->mass;
		r_cm = (l->mass * l->r_cm + r->mass * r->r_cm) / mass;
	}
	else {
		// should never run
		std::cout << "Check MassTree.ccp: Tree::computeMassMoments method." << std::endl;
		if (l == nullptr) {
			mass = r->mass;
			r_cm = r->r_cm;
		}
		else {
			mass = l->mass;
			r_cm = l->r_cm;
		}
	}
}

const int b = 1;
void buildTree(class Tree* T, std::list<class Particle*>& ptrs, std::list<class Particle*>::iterator start, std::list<class Particle*>::iterator stop, std::list<class Particle*>::iterator part,int index, int depth) {
    if (index > 2){
    index=0;
    depth+=1;
  }
  auto left_dist = distance(start, part);
  auto right_dist = distance(part, stop);
  //std::cout<<left_dist<<std::endl;
  if (left_dist+right_dist > 2*b && depth<5) {
    //std::cout<<left_dist<<" , "<<right_dist<<std::endl;
    Tree* left = new Tree(ptrs, start, part,index+1,depth);
    T->l = left;
    if (left_dist > b) {
      std::list<class Particle*>::iterator pl = left->sortParticles(ptrs,start,part,index,depth);
      //std::cout<<"hi";
      buildTree(left, ptrs, start, part,pl,index+1,depth);
    }
    left->computeMassMoments(start, part);
    Tree* right = new Tree(ptrs, part, stop,index+1,depth);
    T->r = right;
    if (right_dist > b) {
      std::list<class Particle*>::iterator pr = right->sortParticles(ptrs,part,stop,index,depth);
      //std::cout<<"hi";
      buildTree(right, ptrs, part, stop, pr,index+1,depth);
    }
    right->computeMassMoments(part, stop);
  }
  //std::cout<<"built"<<depth<<std::endl;
}

void orderParticles(std::vector<class Particle>& ps, std::list<class Particle*>& ptrs) {
	std::vector<class Particle> tmp;
	for (std::list<class Particle*>::iterator p = ptrs.begin(); p != ptrs.end(); p++) tmp.push_back(**p);
	ps.assign(tmp.begin(), tmp.end());
}
