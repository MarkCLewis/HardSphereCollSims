// CollisionTest.cpp
// This is the main file for the rings simulation.

#include<string>
#include<omp.h>

//#define PARALLEL

#include "System.h"
#include "BasicPopulation.h"
#include "CollisionForcing.h"
#include "BoundaryConditions.h"
#include "VariableGridCollisionHash.h"
#include "FixedGridCollisionHash.h"
#include "Distributions.h"

int main(int argc,char **argv) {

/***** Forcing Setup ********/
	typedef FixedGridCollisionHash Hash;
//	typedef VariableGridCollisionHash Hash;
//	typedef CollisionForcing<Hash> CollForcing;
	typedef CollisionForcing<Hash> Forcing;
//	typedef DoubleForce<MoonForcing,CollForcing> Forcing;

//	CollForcing collForce;
	Forcing force;
//	MoonForcing moon(moonMass,moonE,moonPhi);
//	Forcing force(moon,collForce);


/***** Boundary Setup ********/
//	typedef SlidingBrick Boundary;
	typedef PeriodicWithPhiShift Boundary;
	Boundary bc(minx,maxx,miny,maxy);

//	typedef SingleOrbitAzimuthal Boundary;
//	Boundary bc(minx,maxx,miny);

//	typedef ParallelPeriodicWithPhiShift<CollForcing> Boundary;
//	Boundary bc(minx,maxx,miny,maxy,pc,collForce);

//	typedef ParallelSingleOrbitAzimuthal<CollForcing> Boundary;
//	Boundary bc(minx,maxx,miny,pc,collForce);


/***** Output Setup ********/
	typedef VelocityEllipsoidOutput Output;
	Output output;

//	typedef BinaryDumpOutput Output;
//	Output output(outputInterval);

//	typedef ParallelBinaryDumpOutput Output;
//	Output output(outputInterval,pc);

//	typedef DoubleOutput<???,???> Output;


/***** Population Setup ********/
	typedef BasicPopulation<Boundary,LinearFinder> Pop;
	Pop pop(bc);


/***** Particle Distribution Setup ********/
	RadiusDistrib rd(singleRadius);
//	RandomSquareEIOrbits distrib(numBodies,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd);
//	RandomSquareEIOrbitsWithCheck<Boundary> distrib(numBodies,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd,bc);
	TauRandomSquareEIOrbitsWithCheck<Boundary> distrib(tau,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd,bc);
//	TauRandomGaussianEIOrbitsWithCheck<Boundary> distrib(tau,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd,bc);
	pop.randomDistribution(distrib);

	printf("In main %e %e\n",pop.getx(0),pop.gety(0));

/***** Do Simulation ********/
	System<Pop,Forcing> sys(pop,force);

	printf("Start Simulation with %d particles.\n",pop.getNumBodies());
	for(int i=0; pop.getMiny()>stopAzimuth; ++i) {
		printf("Step %d\n",i);
		sys.advance();
	}

	return 0;
}
