// RingsSim.cpp
// This is the main file for the rings simulation.

#include<string>

//#define PARALLEL

#include "CollisionFinders.h"
#include "System.h"
#include "GCPopulation.h"
//#include "TreeCollisionForcing.h"
#include "CollisionForcing.h"
#include "BoundaryConditions.h"
#include "VariableGridCollisionHash.h"
#include "FixedGridCollisionHash.h"
#include "Distributions.h"
#include "VelocityEllipsoidOutput.h"
#include "BinaryDumpOutput.h"
#include "ParallelBinaryDumpOutput.h"
#include "MoonForcing.h"
#include "DoubleForce.h"
#include "ParticleMoonForcing.h"
//#include "PSKDTree.h"

int main(int argc,char **argv) {
//	double minx=0.0060,maxx=0.0068,miny=0.04,maxy=0.1022;
	double minx=0.0066,maxx=0.0068,miny=0.04,maxy=0.0402;
//	double minx=0.0065,maxx=0.006501,miny=0.04,maxy=0.04001;
	double moonMass=2.5e-10;
	double moonE=0.0024;
	double moonPhi=0.0;
	double singleRadius=1e-7;
	int numBodies=10000;
	double tau=0.004;
	double eValue=1e-9;
	double iValue=2e-8;
	double stopAzimuth=-100.2;
	int outputInterval=100;
	int step;
	
	if(argc<2) {
		printf("You have to enter a step number that you are recovering from.");
		return -1;
	}
	step=atoi(argv[1]);
	

/***** Forcing Setup ********/
	typedef FixedGridCollisionHash Hash;
//	typedef VariableGridCollisionHash Hash;
	typedef CollisionForcing<Hash> CollForcing;
	typedef DoubleForce<MoonForcing,CollForcing> Forcing1;
	typedef DoubleForce<ParticleMoonForcing,Forcing1> Forcing;
//	typedef DoubleForce<ParticleMoonForcing,CollForcing> Forcing;

	CollForcing collForce;
//	Forcing force;
	MoonForcing moon(moonMass,moonE,moonPhi);
	Forcing1 force1(moon,collForce);
	ParticleMoonForcing pmf(1e-7,ParticleMoonForcing::calcSaturnianDensity(1.4e5,0.6));
	Forcing force(pmf,force1);
//	Forcing force(pmf,collForce);


/***** Boundary Setup ********/
//	typedef SlidingBrick Boundary;
//	typedef PeriodicWithPhiShift Boundary;
//	Boundary bc(minx,maxx,miny,maxy);

	typedef SingleOrbitAzimuthal Boundary;
	miny-=2.0*GCCoords::A0*minx*5e-3*(step+1);
	Boundary bc(minx,maxx,miny);

//	typedef ParallelPeriodicWithPhiShift<CollForcing> Boundary;
//	Boundary bc(minx,maxx,miny,maxy,pc,collForce);

//	typedef ParallelSingleOrbitAzimuthal<CollForcing> Boundary;
//	Boundary bc(minx,maxx,miny,pc,collForce);


/***** Output Setup ********/
//	typedef VelocityEllipsoidOutput Output;
//	Output output;

	typedef BinaryDumpOutput<CartCoords> Output;
	Output output(outputInterval,step);

//	typedef ParallelBinaryDumpOutput Output;
//	Output output(outputInterval,pc);

//	typedef DoubleOutput<???,???> Output;


/***** Population Setup ********/
	typedef GCPopulation<Boundary,FullNewtonFinder,Output> Pop;
	Pop pop(bc,output);


/***** Particle Distribution Setup ********/
	RadiusDistrib rd(singleRadius);
//	RandomSquareEIOrbits distrib(numBodies,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd);
//	RandomSquareEIOrbitsWithCheck<Boundary> distrib(numBodies,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd,bc);
//	TauRandomSquareEIOrbitsWithCheck<Boundary> distrib(tau,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd,bc);
//	TauRandomGaussianEIOrbitsWithCheck<Boundary> distrib(tau,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd,bc);
	FileRecover distrib(step,pop.getTimeStep(),bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY());
	pop.randomDistribution(distrib);
	
//	pop.addSingleParticleGC((bc.getMaxX()+bc.getMinX())*0.5,(bc.getMaxY()+bc.getMinY())*0.5,0.0,0.0,0.0,0.0,1e-6);

	printf("In main %e %e\n",pop.getx(0),pop.gety(0));

/***** Do Simulation ********/
	System<Pop,Forcing> sys(pop,force);

	printf("Start Simulation with %d particles.\n",pop.getNumBodies());
	for(int i=0; pop.getMiny()>stopAzimuth; ++i) {
		step=i;
		printf("Step %d\n",i);
		sys.advance();
	}
	

	return 0;
}
