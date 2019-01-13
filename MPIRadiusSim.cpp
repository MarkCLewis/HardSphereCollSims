// RingsSim.cpp
// This is the main file for the rings simulation.

#include<string>
#include<omp.h>

#define PARALLEL

#define AZIMUTHAL_MIRRORS

#define BUFFER_WIDTH_MULTIPLE 3

#include "CollisionFinders.h"
#include "System.h"
#include "GCPopulation.h"
//#include "TreeCollisionForcing.h"
#include "CollisionForcing.h"
#include "ParallelBoundaryConditions.h"
#include "VariableGridCollisionHash.h"
#include "FixedGridCollisionHash.h"
#include "Distributions.h"
#include "VelocityEllipsoidOutput.h"
#include "BinaryDumpOutput.h"
#include "SpecificFixed2DBinnedOutput.h"
#include "ParallelOutput.h"
#include "Particle2DBinnedOutput.h"
#include "SpecificBinaryDumpOutput.h"
#include "MoonForcing.h"
#include "DoubleForce.h"
#include "SpecificDoubleOutput.h"
#include "ParticleMoonForcing.h"
//#include "PSKDTree.h"
#include "ParallelGravTree.h"

class NoOutput {
public:
        template<class Population>
        void output(Population &pop) {}
};

int main(int argc,char **argv) {
	omp_set_num_threads(1);
	//just to make sure..
//	double moonMass=2.5e-10*0.2; // Prometheus downsized
//	double moonMass=2.5e-10; // Prometheus
//	double moonMass=1.25e-13; // Daphnis
//	double moonE=0;
//	double moonPhi=0.0;
//	int moonCount=30;
	double singleRadius=1e-8/3.0;
//	int numBodies=10000;
	double tau=1.0;
	printf("argc=%d\n",argc);
	if (argc == 3) {
		omp_set_num_threads(atoi(argv[1]));
		tau = atof(argv[2]);
	}
	printf("OpenMP numThreads = %d\n", omp_get_max_threads());
	printf("tau = %lf\n", tau);
	double eValue=1e-9;
	double iValue=2e-8;
//	double stopAzimuth=-3;
	int outputInterval=100/4;
	double timeStep=6.28e-4;

	double miny=-0.0001,maxy=0.0001,minx=-0.0001,maxx=0.0001;

	printf("size and num %e, %e -> %e\n",(maxx-minx),(maxy-miny),(maxy-miny)*(maxx-minx)/(3.14159*singleRadius*singleRadius)*tau);

    ProcessorCommunication pc(argc,argv);
    srand48(pc.getProcessNum());

    printf("In process %d\n",pc.getProcessNum());

//    double start=MPI_Wtime();	

/***** Forcing Setup ********/
	typedef FixedGridCollisionHash Hash;
//	typedef VariableGridCollisionHash Hash;
	typedef CollisionForcing<Hash> Forcing;
//	typedef DoubleForce<MoonForcing,CollForcing> Forcing1;
		
	typedef SpecificFixedPeriodic SpecificBC;
	SpecificBC sbc(false);
	
	typedef ParallelBoundary<SpecificBC,Forcing> Boundary;
	
//	typedef DoubleForce<KDGravTree<Boundary>,Forcing1> Forcing;
//	typedef DoubleForce<ParticleMoonForcing,CollForcing> Forcing;

	Forcing force;
	
	Boundary bc(sbc,minx,maxx,miny,maxy,pc,force);
	
//	Forcing force;
//	MoonForcing moon(moonCount,moonX,moonY,moonMass,moonE,moonPhi);
//	ParticleMoonForcing pmf(1e-7,ParticleMoonForcing::calcSaturnianDensity(1.4e5,0.6));
//	KDGravTree<Boundary> tree(0.3,pc,bc,1e-6);
//	Forcing force(tree,force1);
//	Forcing force(pmf,force1);
//	Forcing force(pmf,collForce);


/***** Boundary Setup ********/
//	typedef ParallelPeriodicWithPhiShift<CollForcing> Boundary;
//	Boundary bc(minx,maxx,miny,maxy,pc,collForce);

//	typedef ParallelSingleOrbitAzimuthal<CollForcing> Boundary;
//	Boundary bc(minx,maxx,miny,pc,collForce);

//	typedef SpecificFixedPeriodic SpecificBC;
//	SpecificBC sbc(true);
	
//	typedef ParallelBoundary<SpecificBC,CollForcing> Boundary;
//	Boundary bc(sbc,minx,maxx,miny,maxy,pc,collForce);

/***** Output Setup ********/
//	typedef VelocityEllipsoidOutput Output;
//	Output output;

	typedef SpecificBinaryDumpOutput SpecificOutput1;
	SpecificOutput1 so1(outputInterval*10);

	typedef SpecificFixed2DBinnedOutput<Boundary> SpecificOutput2;
	SpecificOutput2 so2(outputInterval,300,3500,minx-2e-5,maxx,bc);

	typedef SpecificDoubleOutput<SpecificOutput1,SpecificOutput2> SpecificOutput;
	SpecificOutput so(so1,so2);


//	typedef ParallelOutput<SpecificOutput> Output;
//	Output output(pc,so);
	typedef NoOutput Output;
	Output output;

//	typedef ParallelFixed2DBinnedOutput<Boundary> Output2;
//	Output2 output2(1000,pc,200,5000,minx-1e-5,maxx+1e-5,bc);

//	typedef Particle2DBinnedOutput Output3;
//	Output3 output3(50,200,1);

//	typedef ParallelBinaryDumpOutput Output;
//	Output output(outputInterval,pc);

//	typedef DoubleOutput<Output2,Output3> Output4;

//	typedef DoubleOutput<Output1,Output2> Output;
//	Output output(output1,output2);


/***** Population Setup ********/
	typedef GCPopulation<Boundary,FullNewtonFinder,Output> Pop;
	Pop pop(bc,output,timeStep,136530,0.7);


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
	for(int i=0; i<20; ++i) {
		step=i;
		printf("Step %d\n",i);
		sys.advance();
	}
	

	return 0;
}
