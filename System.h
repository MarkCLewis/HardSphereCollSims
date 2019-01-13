// System.h
// This is the main class that I will use.  It will have the main loop in it
// and will get templated on the various components of the simulation.

#ifndef SYSTEM
#define SYSTEM

#include "ParticleIndex.h"

template<class Population, class Forces>
class System {
	public:
		System(Population &p,Forces &f):pop(p),force(f) {}

		void advance() {
			force.applyForce(pop);
			printf("Forcing done\n");
			fflush(stdout);
			pop.endStep();
			printf("Step ended\n");
			fflush(stdout);
		}
	private:
		Population &pop;
		Forces &force;
};

#endif
