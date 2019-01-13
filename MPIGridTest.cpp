// RingsSim.cpp
// This is the main file for the rings simulation.

#include<string>
#include<omp.h>

#define PARALLEL
//#define GRAVITY

//#define AZIMUTHAL_MIRRORS

#define BUFFER_WIDTH_MULTIPLE 3

#include "CollisionFinders.h"
#include "System.h"
#include "GCPopulation.h"
//#include "TreeCollisionForcing.h"
#include "ParallelBoundaryConditions.h"
#include "CollisionForcing.h"
//#include "VariableGridCollisionHash.h"
#include "FixedGridCollisionHash.h"
#include "Distributions.h"
#include "SigmaDistributions.h"
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
//#include "ParallelGravCollTree.h"

class NoOutput {
public:
        template<class Population>
        void output(Population &pop) {}
};

int main(int argc,char **argv) {
	omp_set_num_threads(8);
//	double moonMass=2.5e-10*0.2; // Prometheus downsized
//	double moonMass=2.5e-10; // Prometheus
//	double moonMass=1.25e-13; // Daphnis
//	double moonE=0;
//	double moonPhi=0.0;
//	int moonCount=30;
//	double minx=-0.00001,maxx=0.00001,miny=-0.0003,maxy=0.0003;
	double minx=-0.00001,maxx=0.00001,miny=-0.00003,maxy=0.00003;

	double numBodies = 450000;
	double minRadius=1e-8;//5e-9;
	double maxRadius=5e-8;
	double q=2.8;
	double tau=0.1;
//	double eValue=1e-9;
//	double iValue=2e-8;
	double eValue=1e-9;
	double iValue=2e-8;
	double stopAzimuth=-100.2;
	int outputInterval=100;
	double timeStep=6.28e-3;
	double moonOrbitRadius=130000;
	double particleDensitygPercm3=0.7;
	int timeStepMultiple=4;
	double sigma=40;
	int step=0;

    ProcessorCommunication pc(argc,argv);
    srand48(pc.getProcessNum());

    printf("In process %d\n",pc.getProcessNum());

//    double start=MPI_Wtime();	

/***** Boundary Setup ********/
//	typedef ParallelPeriodicWithPhiShift<CollForcing> Boundary;
//	Boundary bc(minx,maxx,miny,maxy,pc,collForce);

//	typedef ParallelSingleOrbitAzimuthal<CollForcing> Boundary;
//	Boundary bc(minx,maxx,miny,pc,collForce);

//	typedef SpecificFixedPeriodic SpecificBC;
//	SpecificBC sbc(true);
	
//	typedef ParallelBoundary<SpecificBC,CollForcing> Boundary;
//	Boundary bc(sbc,minx,maxx,miny,maxy,pc,collForce);

	typedef SpecificFixedPeriodic SpecificBC;
	SpecificBC sbc(false);
	
/***** Forcing Setup ********/
	typedef FixedGridCollisionHash Hash;
//	typedef VariableGridCollisionHash Hash;
	typedef CollisionForcing<Hash> Forcing;
	typedef ParallelBoundary<SpecificBC,Forcing> Boundary;
//	typedef DoubleForce<MoonForcing,CollForcing> Forcing1;
		
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


/***** Output Setup ********/
//	typedef VelocityEllipsoidOutput Output;
//	Output output;

//	typedef SpecificBinaryDumpOutput SpecificOutput;
//	SpecificOutput so(outputInterval/timeStepMultiple);

//	typedef SpecificFixed2DBinnedOutput<Boundary> SpecificOutput2;
//	SpecificOutput2 so2(outputInterval,300,3500,minx-2e-5,maxx,bc);

//	typedef SpecificDoubleOutput<SpecificOutput1,SpecificOutput2> SpecificOutput;
//	SpecificOutput so(so1,so2);


//	typedef ParallelOutput<SpecificOutput> Output;
//	Output output(pc,so,step);
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
	Pop pop(bc,output,timeStep*timeStepMultiple,moonOrbitRadius,particleDensitygPercm3);


/***** Particle Distribution Setup ********/
	RadiusDistrib rd(minRadius);//,maxRadius,q);
//	RandomSquareEIOrbits distrib(numBodies,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd);
//	RandomSquareEIOrbitsWithCheck<Boundary> distrib(numBodies,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd,bc);
//	SigmaRandomSquareEIOrbitsWithCheck<Boundary> distrib(sigma,particleDensitygPercm3,moonOrbitRadius,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd,bc);
	TauRandomSquareEIOrbitsWithCheck<Boundary> distrib(tau,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd,bc);
//	TauRandomGaussianEIOrbitsWithCheck<Boundary> distrib(tau,bc.getMinX(),bc.getMaxX(),bc.getMinY(),bc.getMaxY(),eValue,iValue,rd,bc);
//	ParallelFileRecover32 distrib(step,1,minx,maxx,miny,maxy,pc);
	pop.randomDistribution(distrib);
	
//	pop.addSingleParticleGC((bc.getMaxX()+bc.getMinX())*0.5,(bc.getMaxY()+bc.getMinY())*0.5,0.0,0.0,0.0,0.0,1e-6);

	printf("In main %e %e\n",pop.getx(0),pop.gety(0));

/***** Do Simulation ********/
	System<Pop,Forcing> sys(pop,force);

	printf("Start Simulation with %d particles.\n",pop.getNumBodies());
	for(int i=0; i<5; ++i) {
		step=i;
		printf("Step %d\n",i);
		sys.advance();
	}
	

	return 0;
}

