// AdhesionForces.h
// This file contains classes for doing different accretion forces.

#ifndef ADHESION_FORCING
#define ADHESION_FORCING

#include<vector>

using std::vector;

struct AdhesionAcc {
	double x,y,z;
};

class NoAdhesionForce {
    public:
        template<class Population,class HashStructure>
        void doAdhesionForce(Population &pop,HashStructure &hash) {}
};

class DoForceAdhesion {
	public:
		DoForceAdhesion() {
			acc.resize(omp_get_max_threads());
		}

		template<class Population>
		void resizeAccs(Population &pop) {
			for(unsigned int t=0; t<acc.size(); ++t) {
				acc[t].resize(pop.getNumBodies());
				for(unsigned int i=0; i<acc[t].size(); ++i) {
					acc[t][i].x = 0.0;
					acc[t][i].y = 0.0;
					acc[t][i].z = 0.0;
				}
			}
		}

		template<class Population>
		void applyAccs(Population &pop) {
			for(unsigned int t=0; t<acc.size(); ++t) {
				for(unsigned int i=0; i<acc[t].size(); ++i) {
					pop.setvx(i,pop.getvx(i)+pop.getTimeStep()*acc[t][i].x);
					pop.setvy(i,pop.getvy(i)+pop.getTimeStep()*acc[t][i].y);
					pop.setvz(i,pop.getvz(i)+pop.getTimeStep()*acc[t][i].z);
				}
			}
			for(int i=0; i<pop.getNumBodies(); ++i) {
				pop.adjustAfterForce(i);
			}
		}
	protected:
		vector<vector<AdhesionAcc> > acc;
};

class SquareForceAdhesion:public DoForceAdhesion {
	public:
		SquareForceAdhesion():forceLength(0),forceMag(0) {}

		/**
		 * This is a fraction beyond the sum of the radii that the force
		 * should act. So if set to 0.1 the attractive force will be felt
		 * when the particles are at 1.1*(r1+r2) or closer.
		 */
		void setForceLength(double fl) {
			forceLength=fl;
		}

		/**
		 * The magnitude is a multiple of the gravity when touching.
		 */
		void setForceMagnitude(double fm) {
			forceMag=fm;
		}

		template<class HashStructure>
		void setupHash(HashStructure &hash) {
			hash.setMinimumGridSize(forceLength);
		}

		template<class Population,class CollisionForce>
		void doAdhesionForce(Population &pop,CollisionForce &cf) {
			resizeAccs(pop);
			cf.runThroughCollisionPairs(pop,*this);
			applyAccs(pop);
		}

		template<class Population>
		void operator()(Population &pop,int p1,int p2) {
			double dx=pop.getx(p2)-pop.getx(p1);
			double dy=pop.gety(p2)-pop.gety(p1);
			double dz=pop.getz(p2)-pop.getz(p1);
			double dist=sqrt(dx*dx+dy*dy+dz*dz);
			if(dist>(pop.getRadius(p1)+pop.getRadius(p2))*(1+forceLength))
				return;
			double fm=calcForceMag(pop,p1,p2)/dist;
			double mag=fm*pop.getMass(p2);
			acc[omp_get_thread_num()][p1].x+=dx*mag;
			acc[omp_get_thread_num()][p1].y+=dy*mag;
			acc[omp_get_thread_num()][p1].z+=dz*mag;
			mag=fm*pop.getMass(p1);
			acc[omp_get_thread_num()][p2].x-=dx*mag;
			acc[omp_get_thread_num()][p2].y-=dy*mag;
			acc[omp_get_thread_num()][p2].z-=dz*mag;
		}

		template<class Population>
		double calcForceMag(Population &pop,int p1,int p2) {
			double dist=pop.getRadius(p1)+pop.getRadius(p2);
			double mag=forceMag/(dist*dist);
			return mag;
		}

	private:
		double forceLength;
		double forceMag;
};

template<class Condition>
class OptionalSquareForceAdhesion:public DoForceAdhesion {
	public:
		OptionalSquareForceAdhesion():forceLength(0),forceMag(0) {}

		/**
		 * This is a fraction beyond the sum of the radii that the force
		 * should act. So if set to 0.1 the attractive force will be felt
		 * when the particles are at 1.1*(r1+r2) or closer.
		 */
		void setForceLength(double fl) {
			forceLength=fl;
		}

		/**
		 * The magnitude is a multiple of the gravity when touching.
		 */
		void setForceMagnitude(double fm) {
			forceMag=fm;
		}

		template<class HashStructure>
		void setupHash(HashStructure &hash) {
			hash.setMinimumGridSize(forceLength);
		}

		template<class Population,class CollisionForce>
		void doAdhesionForce(Population &pop,CollisionForce &cf) {
			resizeAccs(pop);
			cf.runThroughCollisionPairs(pop,*this);
			applyAccs(pop);
		}

		template<class Population>
		void operator()(Population &pop,int p1,int p2) {
			if(!cond(p1,p2)) return;
			double dx=pop.getx(p2)-pop.getx(p1);
			double dy=pop.gety(p2)-pop.gety(p1);
			double dz=pop.getz(p2)-pop.getz(p1);
			double dist=sqrt(dx*dx+dy*dy+dz*dz);
			if(dist>(pop.getRadius(p1)+pop.getRadius(p2))*(1+forceLength))
				return;
			double fm=calcForceMag(pop,p1,p2)/dist;
			double mag=fm*pop.getMass(p2);
			acc[omp_get_thread_num()][p1].x+=dx*mag;
			acc[omp_get_thread_num()][p1].y+=dy*mag;
			acc[omp_get_thread_num()][p1].z+=dz*mag;
			mag=fm*pop.getMass(p1);
			acc[omp_get_thread_num()][p2].x-=dx*mag;
			acc[omp_get_thread_num()][p2].y-=dy*mag;
			acc[omp_get_thread_num()][p2].z-=dz*mag;
		}

		template<class Population>
		double calcForceMag(Population &pop,int p1,int p2) {
			double dist=pop.getRadius(p1)+pop.getRadius(p2);
			double mag=forceMag/(dist*dist);
			return mag;
		}

	private:
		double forceLength;
		double forceMag;
		Condition cond;
};

class SmoothForceAdhesion:public DoForceAdhesion {
	public:
		SmoothForceAdhesion():forceLength(0),forceMag(0) {}

		/**
		 * This is a fraction beyond the sum of the radii that the force
		 * should act. So if set to 0.1 the attractive force will be felt
		 * when the particles are at 1.1*(r1+r2) or closer.
		 */
		void setForceLength(double fl) {
			forceLength=fl;
			b=16/(forceLength*forceLength*forceLength*forceLength);
			c=-32/(forceLength*forceLength*forceLength);
			d=16/(forceLength*forceLength);
		}

		/**
		 * The magnitude is a multiple of the gravity when touching.
		 */
		void setForceMagnitude(double fm) {
			forceMag=fm;
		}

		template<class HashStructure>
		void setupHash(HashStructure &hash) {
			hash.setMinimumGridSize(forceLength);
		}

		template<class Population,class CollisionForce>
		void doAdhesionForce(Population &pop,CollisionForce &cf) {
			resizeAccs(pop);
			cf.runThroughCollisionPairs(pop,*this);
			applyAccs(pop);
		}

		template<class Population>
		void operator()(Population &pop,int p1,int p2) {
			double dx=pop.getx(p2)-pop.getx(p1);
			double dy=pop.gety(p2)-pop.gety(p1);
			double dz=pop.getz(p2)-pop.getz(p1);
			double dist=sqrt(dx*dx+dy*dy+dz*dz);
			if(dist<pop.getRadius(p1)+pop.getRadius(p2) || dist>(pop.getRadius(p1)+pop.getRadius(p2))*(1+forceLength))
				return;
			double fm=calcForceMag(pop,p1,p2,dist)/dist;
			double mag=fm*pop.getMass(p2);
			acc[omp_get_thread_num()][p1].x+=dx*mag;
			acc[omp_get_thread_num()][p1].y+=dy*mag;
			acc[omp_get_thread_num()][p1].z+=dz*mag;
			mag=fm*pop.getMass(p1);
			acc[omp_get_thread_num()][p2].x-=dx*mag;
			acc[omp_get_thread_num()][p2].y-=dy*mag;
			acc[omp_get_thread_num()][p2].z-=dz*mag;
		}

		template<class Population>
		double calcForceMag(Population &pop,int p1,int p2,double pdist) {
			double dist=pop.getRadius(p1)+pop.getRadius(p2);
			pdist = pdist/dist-1;
			double pdsqr=pdist*pdist;
			double mag=forceMag*(b*pdsqr*pdsqr+c*pdist*pdsqr+d*pdsqr)/(dist*dist);
			return mag;
		}

	private:
		double forceLength;
		double forceMag;
		double b,c,d;
};

template<class Cond>
class ConditionalSmoothForceAdhesion:public DoForceAdhesion {
	public:
		ConditionalSmoothForceAdhesion():forceLength(0),forceMag(0) {}

		/**
		 * This is a fraction beyond the sum of the radii that the force
		 * should act. So if set to 0.1 the attractive force will be felt
		 * when the particles are at 1.1*(r1+r2) or closer.
		 */
		void setForceLength(double fl) {
			forceLength=fl;
			b=16/(forceLength*forceLength*forceLength*forceLength);
			c=-32/(forceLength*forceLength*forceLength);
			d=16/(forceLength*forceLength);
		}

		/**
		 * The magnitude is a multiple of the gravity when touching.
		 */
		void setForceMagnitude(double fm) {
			forceMag=fm;
		}

		template<class HashStructure>
		void setupHash(HashStructure &hash) {
			hash.setMinimumGridSize(forceLength);
		}

		template<class Population,class CollisionForce>
		void doAdhesionForce(Population &pop,CollisionForce &cf) {
			resizeAccs(pop);
			cf.runThroughCollisionPairs(pop,*this);
			applyAccs(pop);
		}

		template<class Population>
		void operator()(Population &pop,int p1,int p2) {
			if(!Cond::applyForce(pop,p1,p2)) return;
			double dx=pop.getx(p2)-pop.getx(p1);
			double dy=pop.gety(p2)-pop.gety(p1);
			double dz=pop.getz(p2)-pop.getz(p1);
			double dist=sqrt(dx*dx+dy*dy+dz*dz);
			if(dist<pop.getRadius(p1)+pop.getRadius(p2) || dist>(pop.getRadius(p1)+pop.getRadius(p2))*(1+forceLength))
				return;
			double fm=calcForceMag(pop,p1,p2,dist)/dist;
			double mag=fm*pop.getMass(p2);
			acc[omp_get_thread_num()][p1].x+=dx*mag;
			acc[omp_get_thread_num()][p1].y+=dy*mag;
			acc[omp_get_thread_num()][p1].z+=dz*mag;
			mag=fm*pop.getMass(p1);
			acc[omp_get_thread_num()][p2].x-=dx*mag;
			acc[omp_get_thread_num()][p2].y-=dy*mag;
			acc[omp_get_thread_num()][p2].z-=dz*mag;
		}

		template<class Population>
		double calcForceMag(Population &pop,int p1,int p2,double pdist) {
			double dist=pop.getRadius(p1)+pop.getRadius(p2);
			pdist = pdist/dist-1;
			double pdsqr=pdist*pdist;
			double mag=forceMag*(b*pdsqr*pdsqr+c*pdist*pdsqr+d*pdsqr)/(dist*dist);
			return mag;
		}

	private:
		double forceLength;
		double forceMag;
		double b,c,d;
};

#endif
