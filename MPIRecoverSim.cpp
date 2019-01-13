// RingsSim.cpp
// This is the main file for the rings simulation.

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

int main(int argc,char **argv) {
	double minx=0.0065,maxx=0.0067,miny=0.04,maxy=0.1022;
//	double minx=0.0065,maxx=0.0067,miny=0.04,maxy=0.0402;

	ProcessorCommunication pc(argc,argv);
	srand48(pc.getProcessNum());

	printf("In process %d\n",pc.getProcessNum());

	double start=MPI_Wtime();

	typedef CollisionForcing<FixedGridCollisionHash> CollForcing;
	typedef DoubleForce<MoonForcing,CollForcing> Forcing;
	typedef ParallelSingleOrbitAzimuthal<CollForcing> Boundary;
//	typedef ParallelPeriodicWithPhiShift<CollForcing> Boundary;
//	typedef PeriodicWithPhiShift Boundary;
	typedef ParallelBinaryDumpOutput Output;
	typedef GCPopulation<Boundary,FullNewtonFinder,Output> Pop;

	CollForcing collForce;
	MoonForcing moon(5e-10,0.0024,0.0);
	Forcing force(moon,collForce);

	Output pbdo(500,195500,pc);

//	Boundary bc(minx,maxx,miny,maxy,pc,collForce);
	miny-=2.0*GCCoords::A0*minx*1e-2*195501;
	Boundary bc(minx,maxx,miny,pc,collForce);
	Pop pop(bc,pbdo);

	RadiusDistrib rd(1e-7);
//	RandomSquareEIOrbits distrib(10000000/pc.getNumProcesses(),bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),1e-9,1e-7,rd);
//	TauRandomSquareEIOrbitsWithCheck<Boundary> distrib(0.005,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),1e-9,1e-7,rd,bc);
	FileRecover distrib(195500,pop.getTimeStep(),bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY());
	pop.randomDistribution(distrib);

	System<Pop,Forcing> sys(pop,force);

	for(int i=0; pop.getMiny()>-128.0; i++) {
		if(pc.getProcessNum()==0) printf("Step %d\n",i);
		sys.advance();
		if(pc.getProcessNum()==0 && i%500==10) {
			system("rcp *.bin hyperion:/hdb1/Rings/ParallelSim");
			system("rm *.bin");
		}
	}

	double end=MPI_Wtime();

	printf("Run took %e seconds\n",end-start);

	return 0;
}
