/**
 * This file is intended to provide tests for various parts of the simulation
 * code. The way I write my tests they will be silent if they pass and print
 * error messages if they don't.
 */

#include <iostream>
#include <omp.h>
#include <cmath>

#include "ParticleIndex.h"
#include "BoundaryConditions.h"
#include "GCPopulation.h"
#include "MassFunctions.h"
#include "CollisionFinders.h"
#include "TextOutput.h"

double moonOrbitRadius = 130000;
double particleDensitygPercm3 = 0.5;
double timeStep = 1e-3;

template<class Pop>
void printPartV(Pop &pop, int i) {
	ParticleIndex pi{i};
	std::cout << i << ") " << pop.getvx(pi) << " " << pop.getvy(pi) << " " <<
		pop.getvz(pi) << "\n";
}

void assertEquals(const char *const msg, double expected, double actual, double tol = 1e-8) {
	if (abs(expected - actual) > tol) {
		std::cout << msg << "\n";
		std::cout << "Expected " << expected << " but got " << actual << "\n";
	}
}

void headOnCollision(double v1, double r1, double v2, double r2, double vp1, double vp2) {
	typedef OpenBounds Boundary;
	Boundary bc;
	typedef GCPopulation<Boundary, FullLinearFinder<GCCoords>, TextOutput, GCCoords, En_0p5_Et_0p5> Pop;
	StandardMass massFunc(moonOrbitRadius, particleDensitygPercm3);
	TextOutput output(100);
	Pop pop(bc, output, timeStep, moonOrbitRadius, massFunc, particleDensitygPercm3);

	pop.addSingleParticleCart(0, -r1, 0, 0, v1, 0, r1);
	pop.addSingleParticleCart(0, r2, 0, 0, v2, 0, r2);
	ParticleIndex p1{0};
	ParticleIndex p2{1};
	pop.processCollision(p1, p2, 0.0);
	assertEquals("Head on v1:", vp1, pop.getvy(p1));
	assertEquals("Head on v2:", vp2, pop.getvy(p2));
}

int main(int argc, char **argv) {
	headOnCollision(1.0e7, 1.0, -1.0e7, 1.0, -5.0e6, 5.0e6);
}
