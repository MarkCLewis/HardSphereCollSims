/*
 * BasicPolyPopulation.h
 *
 *  Created on: Jun 23, 2010
 *      Author: Cameron Swords
 */

#ifndef BASICPOLYPOPULATION_H_
#define BASICPOLYPOPULATION_H_


#include <vector>
#include "Coordinates.h"

template<class BoundaryCondition,class OutputMethod>
class BasicPolyPopulation {

public:

	// p1, p2 -> indicies of particles
	// returns t -> time since last timestep
	double collisionTime(int p1,int p2) { }

	// p1, p2 -> indicies of particles
	// t -> time since last timestep
	void processCollision(int p1,int p2,double t) { }

private:
	int numBodies;
	int numReal;
	std::vector<BasicCartCoords> cart;
	std::vector<double> radius;
	std::vector<double> time;
	BoundaryCondition &bounds;
	OutputMethod &out;

	double dt;
	double maxRadius;

	long checkCnt;

};


#endif /* BASICPOLYPOPULATION_H_ */
