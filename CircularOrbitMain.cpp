// CircularOrbitMain.cpp
// This is the main file for testing circular orbits in a solar-system type of
// setup. Can be used for speed tests and limited validity of gravity trees.

#include<string>
#include<omp.h>
#include<sys/time.h>

#define GRAVITY

#include "System.h"
#include "BasicPopulation.h"
#include "TreeCollisionForcing.h"
// #include "CollisionForcing.h"
#include "BoundaryConditions.h"
#include "VariableGridCollisionHash.h"
#include "FixedGridCollisionHash.h"
#include "Distributions.h"
#include "BinaryDumpOutput.h"
#include "DoubleForce.h"
#include "DoubleOutput.h"
#include "GravCollTree.h"

class NoOutput {
public:
        template<class Population>
        void output(Population &pop) {}
};

int main(int argc,char **argv) {
	int numBodies = (argc<2)?100:atoi(argv[1]);
	double stopTime = 6.28;
	int outputInterval=20;
	double dt=(argc<3)?1e-3:(atof(argv[2]));
  double velcmPerS = 1.0; // There shouldn't be collisions in this test, so this shouldn't matter.

/***** Boundary Setup ********/
	typedef OpenBounds Boundary;
	Boundary bc;

/***** Forcing Setup ********/
  // typedef TreeCollisionForcing<GravCollTree<Boundary> > CollForcing;
  // typedef DoubleForce<GravCollTree<Boundary>,CollForcing> Forcing;

	typedef GravCollTree<Boundary> Forcing;
  GravCollTree<Boundary> force(0.3,bc);
  // CollForcing collForce(gt);
  // Forcing force(gt,collForce);


/***** Output Setup ********/
//	typedef BinaryDumpOutput<BasicCartCoords> Output;
//	Output output(outputInterval);
	typedef NoOutput Output;
	Output output;
	
/***** Population Setup ********/
  StandardMass mf(1.5e8, 1.41, 2e30);
	typedef BasicPopulation<Boundary,Output> Pop;
	Pop pop(bc, output, dt, velcmPerS, mf);

/***** Particle Distribution Setup ********/
  pop.addSingleParticle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/215.032);
	for (int i = 0; i < numBodies; ++i) {
		double d = 0.1 + i * 5.0 / numBodies;
		double v = sqrt(1.0 / d);
		pop.addSingleParticle(d, 0.0, 0.0, 0.0, v, 0.0, 1e-14);
		pop.addSingleParticle(-d, 0.0, 0.0, 0.0, -v, 0.0, 1e-14);
	}
	
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
