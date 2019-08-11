#ifndef MASS_FUNCTIONS
#define MASS_FUNCTIONS

#include<algorithm>

/**
 * This file contains functions that calculate mass from density and radius.
 * Passing this in as a function allows me to not have conditionals in the
 * most common case and to vary the mass calculation more generally without
 * messing up the case in the population.
 */

class StandardMass {
		double rhoFact;
	public:
		explicit StandardMass(double r0, double densitygPcm3=0.5, double centralMass=5.68e26) {
			r0 *= 1000; // km to m conversion
			rhoFact=1.33333*3.14159*r0*r0*r0*1e3*densitygPcm3/centralMass;
		}
		double operator()(ParticleIndex pi, double r) const {
			return r*r*r*rhoFact;
		}
};

class AlternateRhoIndexMass {
		int altIndex;
		double rhoFact, rhoFact2;
	public:
		AlternateRhoIndexMass(double r0, int alternateRhoIndex, double densitygPcm3, double alternateDensitygPcm3, double centralMass=5.68e26): altIndex(alternateRhoIndex) {
			r0 *= 1000; // km to m conversion
			rhoFact=1.33333*3.14159*r0*r0*r0*1e3*densitygPcm3/centralMass;
			rhoFact2=1.33333*3.14159*r0*r0*r0*1e3*alternateDensitygPcm3/centralMass;
		}
		double operator()(ParticleIndex pi, double r) const {
			if(r>=altIndex) return r*r*r*rhoFact2;
			else return r*r*r*rhoFact;
		}
};

class AlternateRhoRadiusMass {
		int altRadius;
		double rhoFact, rhoFact2;
	public:
		AlternateRhoRadiusMass(double r0, int alternateRhoRadius, double densitygPcm3, double alternateDensitygPcm3, double centralMass=5.68e26): altRadius(alternateRhoRadius) {
			r0 *= 1000; // km to m conversion
			rhoFact=1.33333*3.14159*r0*r0*r0*1e3*densitygPcm3/centralMass;
			rhoFact2=1.33333*3.14159*r0*r0*r0*1e3*alternateDensitygPcm3/centralMass;
		}
		double operator()(ParticleIndex pi, double r) const {
			if(r>=altRadius) return r*r*r*rhoFact2;
			else return r*r*r*rhoFact;
		}
};

template<typename MF>
class MassRamp {
		MF massFunc;
		int &step;
		int fullStep;
	public:
		MassRamp(MF mf, int &s, int full): massFunc(mf), step(s), fullStep(full) {}
		double operator()(ParticleIndex pi, double r) const {
			double m = massFunc(pi, r);
			if(step >= fullStep) return m;
			else return std::max(step, 1)*m/fullStep;
		}
};

#endif
