// RingsSim.cpp
// This is the main file for the rings simulation.

#include <vector>
using namespace std;

#define PARALLEL

#include "CollisionFinders.h"
#include "System.h"
#include "GCPopulation.h"
#include "CollisionForcing.h"
#include "BoundaryConditions.h"
#include "FixedGridCollisionHash.h"
#include "Distributions.h"
#include "ParallelBinaryDumpOutput.h"
#include "MoonForcing.h"
#include "DoubleForce.h"
#include "ParallelFixed2DBinnedOutput.h"
#include "GravTree.h"

int main(int argc,char **argv) {
//	double minx=0.0061,maxx=0.0069,miny=0.04,maxy=0.1022;
//	double minx=0.0065,maxx=0.0067,miny=0.04,maxy=0.1022;
	double minx=0.00012,maxx=0.00042,miny=-0.02,maxy=0.0022;


	ProcessorCommunication pc(argc,argv);
	srand48(pc.getProcessNum());

	printf("In process %d\n",pc.getProcessNum());

	double start=MPI_Wtime();

	typedef CollisionForcing<FixedGridCollisionHash> CollForcing;
	typedef DoubleForce<MoonForcing,CollForcing> Forcing;
//	typedef ParallelSingleOrbitAzimuthal<CollForcing> Boundary;
	typedef ParallelFixedPeriodic<CollForcing> Boundary;
//	typedef ParallelPeriodicWithPhiShift<CollForcing> Boundary;
//	typedef PeriodicWithPhiShift Boundary;
//	typedef ParallelBinaryDumpOutput Output;
	typedef ParallelFixed2DBinnedOutput Output;
	typedef GCPopulation<Boundary,FullNewtonFinder,Output> Pop;

	CollForcing collForce;
	MoonForcing moon(1.25e-13,2e-5,0.0);
//	Forcing force(moon,collForce);
	
	KDGravTree<Pop>(.2) gravForce;
//	Output pbdo(100,pc);
	Output pbdo(100,pc,400,4000,minx,maxx,miny,maxy);

	Boundary bc(minx,maxx,miny,maxy,pc,collForce,true);
//	Boundary bc(minx,maxx,miny,pc,collForce);
	Pop pop(bc,pbdo);

	RadiusDistrib rd(1e-7);
	RandomSquareEIOrbits distrib(10000000/pc.getNumProcesses(),bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),1e-9,1e-7,rd);
//	TauRandomSquareEIOrbitsWithCheck<Boundary> distrib(0.1,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),1e-9,1e-7,rd,bc);
	pop.randomDistribution(distrib);
        
	//System<Pop,Forcing> sys(pop,force);
	System<Pop,KDGravTree> sys(pop,gravForce);
	for(int i=0; pop.getMiny()>-128.0; i++) {
		if(pc.getProcessNum()==0) printf("Step %d\n",i);
		sys.advance();
//		if(pc.getProcessNum()==0 && i%100==5) {
//			system("./Process.sh &");
//		}
	}

	double end=MPI_Wtime();

	printf("Run took %e seconds\n",end-start);

	return 0;
}
