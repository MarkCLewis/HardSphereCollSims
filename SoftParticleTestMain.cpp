// SoftParticleTestMain.cpp
// This is the main file for a small test of short range forces. The main goal of this is to validate
// that the soft-sphere interactions are behaving reasonably for the range of inputs that we expect.

#include<string>
#include<omp.h>
#include<sys/time.h>

#define GRAVITY

#include "System.h"
#include "BasicPopulation.h"
#include "BoundaryConditions.h"
#include "GravCollTree.h"
#include "ShortRangeForces.h"

class NoOutput {
  public:
    template<class Population>
    void output(Population &pop) {}
};

int main(int argc,char **argv) {
	int numVels = 6;
  double vels[] = {1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6};
  int numSizes = 5;
  double sizes[] = {1e-9, 3e-9, 1e-8, 3e-8, 1e-7};
	double dt=(argc<2) ? 1e-4 : (atof(argv[1]));
  double velcmPerS = 1.0; // There shouldn't be collisions in this test, so this shouldn't matter.

/***** Boundary Setup ********/
	typedef OpenBounds Boundary;
	Boundary bc;

/***** Forcing Setup ********/
	typedef GravCollTree<Boundary, ShortRangeGravitySpringDiscontinuous> Forcing;
  ShortRangeGravitySpringDiscontinuous srf;
  Forcing force(0.3,bc, srf);

/***** Output Setup ********/
	typedef NoOutput Output;
	Output output;
	
/***** Population Setup ********/
  StandardMass mf(130000.0);
	typedef BasicPopulation<Boundary,Output> Pop;
	Pop pop(bc, output, dt, velcmPerS, mf);

/***** Particle Distribution Setup ********/
  for (int j = 0; j < numSizes; ++j) {
	  for (int i = 0; i < numVels; ++i) {
      double y = j * numVels + i;
      double r = sizes[j];
      double v = vels[i];
      double x = r + v * 5 * dt;
      pop.addSingleParticle(x, y, 0.0, -v, 0.0, 0.0, r);
      pop.addSingleParticle(-x, y, 0.0, v, 0.0, 0.0, r);
    }
	}
  int numBodies = pop.getNumBodies();
  printf("Particle locations:\n"); 
	for (int i = 0; i < numBodies; ++i) {
		ParticleIndex pi{i};
		printf("%d, %e %e %e, %e %e %e\n", i, pop.getx(pi), pop.gety(pi), pop.getz(pi), pop.getvx(pi), pop.getvy(pi), pop.getvz(pi));
	}
	
/***** Do Simulation ********/
	System<Pop,Forcing> sys(pop,force);

	double start = omp_get_wtime();

	printf("Start Simulation with %d particles.\n",pop.getNumBodies());
	for(int i=0; i < 1000; ++i) {
		printf("Step %d\n",i);
		sys.advance();
    printf("vx[0] = %e\n", pop.getvx(ParticleIndex{0}));
	}

	double end = omp_get_wtime();

	printf("Finished in %g seconds.\n", end-start);
	printf("Particle locations:\n"); 
	for (int i = 0; i < numBodies; ++i) {
		ParticleIndex pi{i};
		printf("%d, %e %e %e, %e %e %e\n", i, pop.getx(pi), pop.gety(pi), pop.getz(pi), pop.getvx(pi), pop.getvy(pi), pop.getvz(pi));
	}

	return 0;
}
