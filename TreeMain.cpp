// RingsSim.cpp
// This is the main file for the rings simulation.

#include<string>
#include<omp.h>

//#define PARALLEL

#include "CollisionFinders.h"
#include "System.h"
#include "GCPopulation.h"
#include "TreeCollisionForcing.h"
//#include "CollisionForcing.h"
#include "BoundaryConditions.h"
//#include "VariableGridCollisionHash.h"
//#include "FixedGridCollisionHash.h"
#include "Distributions.h"
#include "VelocityEllipsoidOutput.h"
#include "BinaryDumpOutput.h"
#include "Fixed2DBinnedOutput.h"
#include "ParallelBinaryDumpOutput.h"
#include "MoonForcing.h"
#include "DoubleForce.h"
#include "ParticleMoonForcing.h"
//#include "PSKDTree.h"
#include "KDTree2.h"

class NoOutput {
public:
	template<class Population>
        void output(Population &pop) {}
};

int main(int argc,char **argv) {
//	double minx=0.0060,maxx=0.0068,miny=0.04,maxy=0.1022;
//	double minx=0.0066,maxx=0.0068,miny=0.04,maxy=0.0402;
//	double minx=0.0065,maxx=0.006501,miny=0.04,maxy=0.04001;
//	double minx=0.00012,maxx=0.00022,miny=-0.02,maxy=0.0022;
	double minx=0.000,maxx=0.00001,miny=-0.00001,maxy=0.00001;

//	double moonMass=2.5e-10;
	double moonMass=1.25e-13;
	double moonE=2e-5;
	double moonPhi=0.0;
	double singleRadius=1e-8;
	int numBodies=10000;
	double tau=1.0;
	double eValue=1e-9;
	double iValue=2e-8;
	double stopAzimuth=-100.2;
	int outputInterval=100;
	double timeStep=6.28e-3;
	

/***** Forcing Setup ********/
//	typedef FixedGridCollisionHash Hash;
//	typedef VariableGridCollisionHash Hash;
//	typedef CollisionForcing<Hash> Forcing;
//	typedef DoubleForce<MoonForcing,CollForcing> Forcing;
//	typedef DoubleForce<ParticleMoonForcing,Forcing1> Forcing;
//	typedef DoubleForce<ParticleMoonForcing,CollForcing> Forcing;
	typedef KDTree Tree;
	typedef TreeCollisionForcing<Tree> Forcing;

	Tree tree;
	Forcing force(tree);
//	Forcing force;
//	MoonForcing moon(moonMass,moonE,moonPhi);
//	Forcing force(moon,collForce);
//	ParticleMoonForcing pmf(1e-7,ParticleMoonForcing::calcSaturnianDensity(1.4e5,0.6));
//	Forcing force(pmf,force1);
//	Forcing force(pmf,collForce);


/***** Boundary Setup ********/
//	typedef SlidingBrick Boundary;
//	typedef PeriodicWithPhiShift Boundary;
	typedef FixedPeriodic Boundary;
	Boundary bc(minx,maxx,miny,maxy);

//	typedef SingleOrbitAzimuthal Boundary;
//	Boundary bc(minx,maxx,miny);

//	typedef ParallelPeriodicWithPhiShift<CollForcing> Boundary;
//	Boundary bc(minx,maxx,miny,maxy,pc,collForce);

//	typedef ParallelSingleOrbitAzimuthal<CollForcing> Boundary;
//	Boundary bc(minx,maxx,miny,pc,collForce);


/***** Output Setup ********/
	//typedef VelocityEllipsoidOutput Output;
	typedef NoOutput Output;
	Output output;

//	typedef BinaryDumpOutput<CartCoords> Output;
//	Output output(outputInterval);

//	typedef Fixed2DBinnedOutput Output;
//	Output output(outputInterval,300,4000,minx-2e-5,maxx);

//	typedef ParallelBinaryDumpOutput Output;
//	Output output(outputInterval,pc);

//	typedef DoubleOutput<???,???> Output;


/***** Population Setup ********/
	typedef GCPopulation<Boundary,FullNewtonFinder,Output> Pop;
	Pop pop(bc,output,timeStep,100000);


/***** Particle Distribution Setup ********/
	RadiusDistrib rd(singleRadius);
//	RandomSquareEIOrbits distrib(numBodies,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd);
//	RandomSquareEIOrbitsWithCheck<Boundary> distrib(numBodies,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd,bc);
	TauRandomSquareEIOrbitsWithCheck<Boundary> distrib(tau,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd,bc);
//	TauRandomGaussianEIOrbitsWithCheck<Boundary> distrib(tau,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd,bc);
	pop.randomDistribution(distrib);
	
//	pop.addSingleParticleGC((bc.getMaxX()+bc.getMinX())*0.5,(bc.getMaxY()+bc.getMinY())*0.5,0.0,0.0,0.0,0.0,1e-6);

	printf("In main %e %e\n",pop.getx(0),pop.gety(0));

/***** Do Simulation ********/
	System<Pop,Forcing> sys(pop,force);

	printf("Start Simulation with %d particles.\n",pop.getNumBodies());
	for(int i=0; i < 10 && pop.getMiny()>stopAzimuth; ++i) {
		step=i;
		printf("Step %d\n",i);
		sys.advance();
	}

	return 0;
}
