// Distributions.h
// This file contains some classes that can produce different distributions
// of particles to be used with the populations.

#ifndef SIGMA_DISTRIBUTIONS
#define SIGMA_DISTRIBUTIONS

#include <sys/types.h>
#include <unistd.h>

#include "GCPopulation.h"

// Saturn
//#define PLANET_MASS_IN_GRAMS 5.6846e29

// Sun
#define PLANET_MASS_IN_GRAMS 2e33 

using std::vector;

template <class PositionCheck,class GCType>
class SigmaRandomSquareEIOrbitsWithCheck {
	public:
		SigmaRandomSquareEIOrbitsWithCheck(double sig,double p,double R_0,double minX,double maxX,double minY,double maxY,double maxE,double maxI,RadiusDistrib &rad,PositionCheck &pc):sigma(sig),rho(p),minx(minX),miny(minY),maxx(maxX),maxy(maxY),maxe(maxE),maxi(maxI),rd(rad),posCheck(pc) {
			curMass=0.0;
			double totArea=pc.getArea();
			R_0*=1e5;
			totMass=totArea*sigma*(R_0*R_0)/PLANET_MASS_IN_GRAMS;
			rho*=(R_0*R_0*R_0)/PLANET_MASS_IN_GRAMS;
			printf("area=%e\n",totArea);
		}

		bool usesGC() {
			return true;
		}

		bool moreParticles() {
			return curMass<totMass;
		}

		void setNextParticle(CartCoords &cart,double &rad) {
		}
#ifdef SPIN
		void setNextParticle(CartCoords &cart,double &rad,SpinVector &sp) {
		}
#endif
		void setNextParticle(GCType &gc,double &rad) {
			rad=rd.getRadius();
			do {
				gc.X=minx+drand48()*(maxx-minx);
				gc.Y=miny+drand48()*(maxy-miny);
				curMass+=1.3333333*3.141592654*rad*rad*rad*rho;
			} while(!posCheck.checkIfIn(gc.X,gc.Y));
			gc.e=maxe*drand48();
			gc.i=maxi*drand48();
			gc.phi=6.28*drand48();
			gc.zeta=6.28*drand48();
		}

	private:
		double sigma,curMass,totMass,rho;
		double minx,miny,maxx,maxy;
		double maxe,maxi;
		RadiusDistrib rd;
		PositionCheck &posCheck;
};

template <class PositionCheck,class GCType>
class SigmaRandomGaussianEIOrbitsWithCheck {
	public:
		SigmaRandomGaussianEIOrbitsWithCheck(double sig,double p,double R_0,double minX,double maxX,double minY,double maxY,double maxE,double maxI,RadiusDistrib &rad,PositionCheck &pc):sigma(sig),rho(p),minx(minX),miny(minY),maxx(maxX),maxy(maxY),maxe(maxE),maxi(maxI),rd(rad),posCheck(pc) {
			curMass=0.0;
			double totArea=pc.getArea();
			R_0*=1e5;
			totMass=totArea*sigma*(R_0*R_0)/PLANET_MASS_IN_GRAMS;
			rho*=(R_0*R_0*R_0)/PLANET_MASS_IN_GRAMS;
			printf("area=%e\n",totArea);
		}

		bool usesGC() {
			return true;
		}

		bool moreParticles() {
			return curMass<totMass;
		}

		void setNextParticle(CartCoords &cart,double &rad) {
		}

#ifdef SPIN
		void setNextParticle(CartCoords &cart,double &rad,SpinVector &sp) {
		}
#endif

		void setNextParticle(GCType &gc,double &rad) {
			double midx=(maxx+minx)/2.0;
			double widthx=(maxx-minx)/4.0;
			rad=rd.getRadius();
			do {
				double r1,r2,s;
				do {
					r1=2.0*drand48()-1.0;
					r2=2.0*drand48()-1.0;
					s=r1*r1+r2*r2;
				} while(s>=1.0);
				gc.X=midx+widthx*r1*sqrt(-2.0*log(s)/s);
				gc.Y=miny+drand48()*(maxy-miny);
			} while(!posCheck.checkIfIn(gc.X,gc.Y));
			curMass+=1.3333333*3.141592654*rad*rad*rad*rho;
			gc.e=maxe*drand48();
			gc.i=maxi*drand48();
			gc.phi=6.28*drand48();
			gc.zeta=6.28*drand48();
		}

	private:
		double sigma,curMass,totMass,rho;
		double minx,miny,maxx,maxy;
		double maxe,maxi;
		RadiusDistrib rd;
		PositionCheck &posCheck;
};

#endif
