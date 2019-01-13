// DomingoMain.cpp
// This is the main file for basic granular flow simulations.

#include<string>

//#define PARALLEL

#include "System.h"
#include "BasicPopulation.h"
#include "TreeCollisionForcing.h"
//#include "CollisionForcing.h"
#include "BoundaryConditions.h"
#include "VariableGridCollisionHash.h"
#include "FixedGridCollisionHash.h"
#include "Distributions.h"
#include "BinaryDumpOutput.h"
#include "DoubleForce.h"
#include "DoubleOutput.h"
#include "PSKDTree.h"
#include "KDTree.h"
#include "BinnedOutput.h"

int main(int argc,char **argv) {
	double sizex=1.0,sizey=1.0,sizez=0.1;
	double shear=0.1;
	double initialVel=0.01;
	double singleRadius=1e-2;
	int numBodies=10000;
	int stopSteps=2000;
	int outputInterval=10;
	double dt=1e-1;

        srand48(347236);

	//srand48(3236);
	

/***** Forcing Setup ********/
//	typedef FixedGridCollisionHash Hash;
//	typedef VariableGridCollisionHash Hash;
//	typedef CollisionForcing<Hash> Forcing;
//	typedef quadtree Tree;
	typedef KDTree Tree;
//	typedef PSKDTree Tree;
	typedef TreeCollisionForcing<Tree> Forcing;

	Forcing force;


/***** Boundary Setup ********/
	typedef ShearedPeriodic Boundary;
	Boundary bc(sizex,sizey,sizez,shear);

/***** Output Setup ********/
	typedef BinaryDumpOutput<BasicCartCoords> BinaryDump;
	BinaryDump dump(outputInterval);
	
	typedef BinnedOutput Binned;
	Binned binned(outputInterval,1,10,3,-1.0,1.0);

	typedef DoubleOutput<BinaryDump,BinnedOutput> Output;
	Output output(dump,binned);

/***** Population Setup ********/
	typedef BasicPopulation<Boundary,Output> Pop;
	Pop pop(bc,output,dt);

/***** Particle Distribution Setup ********/
	RadiusDistrib rd(singleRadius);
	CubeDistrib distrib(numBodies,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),bc.getMinZ(),bc.getMaxZ(),initialVel,rd);
	pop.randomDistribution(distrib);
	
/***** Do Simulation ********/
	System<Pop,Forcing> sys(pop,force);

	printf("Start Simulation with %d particles.\n",pop.getNumBodies());
	for(int i=0; i<stopSteps; ++i) {
                printf("Step %d\n",i);
	//	printf("********************Step %d ***********************\n",i);
		sys.advance();
	}
	

	return 0;
}
