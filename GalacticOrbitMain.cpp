// CircularOrbitMain.cpp
// This is the main file for testing a galactic-style setup. 
// Can be used for speed tests and limited validity of gravity trees.

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
	typedef NoOutput Output;
//	Output output(outputInterval);
	Output output;
	
/***** Population Setup ********/
  StandardMass mf(1.5e8, 1.41, 2e30);
	typedef BasicPopulation<Boundary,Output> Pop;
	Pop pop(bc, output, dt, velcmPerS, mf);

/***** Particle Distribution Setup ********/
  pop.addSingleParticle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/2150.032);
	double m = mf(ParticleIndex{0}, pop.getRadius(ParticleIndex{0}));
	double annulusWidth = 5.0 / numBodies;
	double innerRadius = 0.1;
	for (int i = 0; i < numBodies; ++i) {
		double d = innerRadius + i * annulusWidth;
		double theta = drand48() * 6.28;
		double v = sqrt(m / d);  // We want this to be rougly constant, so place equal mass at each annulus.
		pop.addSingleParticle(d * cos(theta), d*sin(theta), 0.0, -v * sin(theta), v * cos(theta), 0.0, 1e-6);
		m += mf(ParticleIndex{i+1}, pop.getRadius(ParticleIndex{i+1}));
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
