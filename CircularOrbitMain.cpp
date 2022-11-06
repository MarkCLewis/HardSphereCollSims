// CircularOrbitMain.cpp
// This is the main file for testing circular orbits in a solar-system type of
// setup. Can be used for speed tests and limited validity of gravity trees.

#include<string>
#include<omp.h>
#include<sys/time.h>

#define GRAVITY

#include "System.h"
#include "BasicPopulation.h"
#include "BoundaryConditions.h"
#include "Distributions.h"
#include "TextOutput.h"
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
	double dt=(argc<3)?1e-3:(atof(argv[2]));
  double velcmPerS = 1.0; // There shouldn't be collisions in this test, so this shouldn't matter.

/***** Boundary Setup ********/
	typedef OpenBounds Boundary;
	Boundary bc;

/***** Forcing Setup ********/
	typedef GravCollTree<Boundary> Forcing;
  GravCollTree<Boundary> force(0.3,bc);

/***** Output Setup ********/
	typedef NoOutput Output;
	Output output;
	// typedef TextOutput Output;
	// Output output(0.1);

	omp_set_num_threads(1);

	
/***** Population Setup ********/
  StandardMass mf(1.5e8, 1.41, 2e30);
	typedef BasicPopulation<Boundary,Output> Pop;
	Pop pop(bc, output, dt, velcmPerS, mf);

/***** Particle Distribution Setup ********/
  pop.addSingleParticle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/215.032);
	for (int i = 0; i < numBodies; ++i) {
		double d = 0.1 + i * 5.0 / numBodies;
		double v = sqrt(1.0 / d);
		double theta = drand48() * 6.28;
		pop.addSingleParticle(d * cos(theta), d * sin(theta), 0.0, -v * sin(theta), v * cos(theta), 0.0, 1e-14);
	}
	printf("Central mass = %f\n", pop.getMass(ParticleIndex{0}));
	printf("Particle locations:\n"); 
	for (int i = 0; i < numBodies; ++i) {
		ParticleIndex pi{i};
		printf("%d %e %e %e\n", i, pop.getx(pi), pop.gety(pi), pop.getz(pi));
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
	// printf("Particle locations:\n");
	// for (int i = 0; i < numBodies; ++i) {
	// 	ParticleIndex pi{i};
	// 	printf("%d %e %e %e\n", i, pop.getx(pi), pop.gety(pi), pop.getz(pi));
	// }

	return 0;
}
