// BasicMain.cpp
// This is the main file for basic granular flow simulations.

#include<string>
#include<omp.h>
#include<sys/time.h>

//#define PARALLEL

#include "System.h"
#include "BasicPopulation.h"
//#include "TreeCollisionForcing.h"
#include "CollisionForcing.h"
#include "BoundaryConditions.h"
#include "VariableGridCollisionHash.h"
#include "FixedGridCollisionHash.h"
#include "Distributions.h"
#include "BinaryDumpOutput.h"
#include "DoubleForce.h"
#include "DoubleOutput.h"
#include "PSKDTree.h"
//#include "KDTree2.h"
#include "BinnedOutput.h"

int main(int argc,char **argv) {
	double sizex=1.0,sizey=1.0,sizez=0.1;
	double shear=0.1;
	double initialVel=0.1;
//	double singleRadius=6e-3;
	double singleRadius=1e-3;
//	int numBodies=88420;
	int numBodies=884200;
	int stopSteps=2001;
//	double stopTime=30.0;
	double stopTime=0.3;
	int outputInterval=200000;
	double dt=(argc<2)?1e-1:(atof(argv[1]));
	

/***** Forcing Setup ********/
	typedef FixedGridCollisionHash Hash;
//	typedef VariableGridCollisionHash Hash;
	typedef CollisionForcing<Hash> Forcing;
//	typedef PSKDTree Tree;
//	typedef KDTree Tree;
//	typedef TreeCollisionForcing<Tree> Forcing;

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
//	TextFileDistrib distrib("bodyIn.txt");
	pop.randomDistribution(distrib);
	
/***** Do Simulation ********/
	System<Pop,Forcing> sys(pop,force);

	double start = omp_get_wtime();

	printf("Start Simulation with %d particles.\n",pop.getNumBodies());
	for(int i=0; i*dt<stopTime; ++i) {
		printf("Step %d\n",i);
		sys.advance();
	}

	double end = omp_get_wtime();

	printf("Finished in %g seconds.\n", end-start);
	

	return 0;
}
