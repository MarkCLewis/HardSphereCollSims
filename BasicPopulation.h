// BasicPopulation.h
// This file describes a class that has a population of particles whose motion
// is along straight lines.

#ifndef BASIC_POPULATION
#define BASIC_POPULATION

#include <vector>
#include <algorithm>
#include "Coordinates.h"
#include "CoefficientsOfRestitution.h"
#include "MassFunctions.h"


template<class BoundaryCondition,class OutputMethod,class CoefRest = En_Bridges_Et_0p5, class MassFunc=StandardMass>
class BasicPopulation {
	public:
		BasicPopulation(BoundaryCondition &bc,OutputMethod &om,double timeStep, double velConv, MassFunc mf):
				bounds(bc),out(om),dt(timeStep),velTocms(velConv),massFunc(mf) {
			numBodies=0;
			numReal=0;
			maxRadius=0.0;
		}

		// For this functions the Distribution class needs to have two basic
		// pieces of functionality.  It needs to have a function that can
		// create a new particle as well as a function that can tell if another
		// particle needs to be generated.  This latter one allows simulations
		// that aim for a certain average optical depth, but it makes dealing
		// with the arrays more difficult.
		template<class Distribution>
		void randomDistribution(Distribution &pd) {
			numBodies=0;
			double vol=0.0;
			while(pd.moreParticles()) {
				numBodies++;
				cart.resize(numBodies);
				radius.resize(numBodies);
				time.resize(numBodies);
				pd.setNextParticle(cart[numBodies-1],radius[numBodies-1]);
				vol+=1.333333*3.14159*radius[numBodies-1]*radius[numBodies-1]*radius[numBodies-1];
				if(numBodies==1 || radius[numBodies-1]>maxRadius)
					maxRadius=radius[numBodies-1];
			}
			printf("Fill factor is %e\n",vol/((getMaxx()-getMinx())*(getMaxy()-getMiny())*(getMaxz()-getMinz())));
			numReal=numBodies;
			checkCnt=0;
		}

		void addSingleParticle(double x,double y,double z,double vx,double vy,double vz,double rad) {
			numBodies++;
			numReal++;
			cart.resize(numBodies);
			radius.resize(numBodies);
			time.resize(numBodies);
			cart[numBodies-1].p[0]=x;
			cart[numBodies-1].p[1]=y;
			cart[numBodies-1].p[2]=z;
			cart[numBodies-1].p[3]=vx;
			cart[numBodies-1].p[4]=vy;
			cart[numBodies-1].p[5]=vz;
			radius[numBodies-1]=rad;
			if(numBodies==1 || radius[numBodies-1]>maxRadius)
				maxRadius=radius[numBodies-1];
		}

		void endStep(void) {
			// Advance all particles to end of step.
			#pragma omp parallel for schedule(static)
			for(int i=0; i<numBodies; ++i) {
				double t=dt-time[i];
				AdvanceParticle(t,ParticleIndex{i});
				time[i]=0.0;
			}

			bounds.apply(*this);
			out.output(*this);
			
			printf("%ld checks\n",checkCnt);
			checkCnt=0;
		}

		double collisionTime(ParticleIndex p1, ParticleIndex p2) {
			checkCnt++;
			double t=(time[p1.i]>time[p2.i])?time[p1.i]:time[p2.i];
			double dt1=t-time[p1.i];
			double dt2=t-time[p2.i];
			double dx=(getx(p2)+dt2*getvx(p2))-(getx(p1)+dt1*getvx(p1));
			double dy=(gety(p2)+dt2*getvy(p2))-(gety(p1)+dt1*getvy(p1));
			double dz=(getz(p2)+dt2*getvz(p2))-(getz(p1)+dt1*getvz(p1));
			double c=dx*dx+dy*dy+dz*dz-(radius[p1.i]+radius[p2.i])*(radius[p1.i]+radius[p2.i]);
			if(c<-1e-4*(radius[p1.i]+radius[p2.i])*(radius[p1.i]+radius[p2.i])) {
//				printf("Overlapping particle %d %d\n",p1,p2);
				if(t>0.0) {
					double dvx=getvx(p2)-getvx(p1);
					double dvy=getvy(p2)-getvy(p1);
					double dvz=getvz(p2)-getvz(p1);
					double vel=sqrt(dvx*dvx+dvy*dvy+dvz*dvz);
					double minRad=(radius[p1.i]<radius[p2.i])?radius[p1.i]:radius[p2.i];
					double time=t+std::min(0.001*minRad/vel,0.5*getTimeStep());
					if(time>getTimeStep()) return -4.0;
					return time;
				} else {
					return 0.0;
				}
			}
			double dvx=getvx(p2)-getvx(p1);
			double dvy=getvy(p2)-getvy(p1);
			double dvz=getvz(p2)-getvz(p1);
			double b=2.0*(dvx*dx+dvy*dy+dvz*dz);
			//if(b>0.0) return -3.0;	// particles heading apart
			double a=dvx*dvx+dvy*dvy+dvz*dvz;
			double rootPart=b*b-4*a*c;
			if(rootPart<0) return -1.0;
			double ret=(-b-sqrt(rootPart))/(2*a)+time[p1.i];
//			printf("%d %d %e %e %e\n",p1,p2,ret,time[p1.i],time[p2.i]);
			if(ret<=time[p1.i] || ret<=time[p2.i]) return -2.0;
			return ret;
		}

		void processCollision(ParticleIndex p1,ParticleIndex p2,double t) {
			double cmpx,cmpy,cmpz,dx,dy,dz,gn,epsilon;
			double gnx,gny,gnz,gtx,gty,gtz;
			double px,py,pz;
			double dmag;
			double mass1=radius[p1.i]*radius[p1.i]*radius[p1.i];
			double mass2=radius[p2.i]*radius[p2.i]*radius[p2.i];
			double tmass=mass1+mass2;
			//double deltax[3];
		
		//	printf("1. Average X=%1.16e\n",0.5*(gc[p1.i].X+gc[p2.i].X));
		
			// Advance the particles to the time of the collision.
			if((t<time[p1.i]) || (t<time[p2.i])) {
				printf("Bad time in collision! %e %e %e %d %d\n",t,time[p1.i],time[p2.i],p1.i,p2.i);
				return;
			}
		
			AdvanceParticle(t-time[p1.i],p1);
			AdvanceParticle(t-time[p2.i],p2);
		
//			printf("Collision between %d and %d at %e\n",p1,p2,t);
		//	printf("Initial vels (%e %e %e)\n   (%e %e %e)\n",cart[p1.i].vx,cart[p1.i].vy,cart[p1.i].vz,cart[p2.i].vx,cart[p2.i].vy,cart[p2.i].vz);
			// Move to the center of mass frame.
			cmpx=(mass2*getvx(p2)+mass1*getvx(p1))/tmass;
			cmpy=(mass2*getvy(p2)+mass1*getvy(p1))/tmass;
			cmpz=(mass2*getvz(p2)+mass1*getvz(p1))/tmass;
			px=(getvx(p1)-cmpx)*mass1;
			py=(getvy(p1)-cmpy)*mass1;
			pz=(getvz(p1)-cmpz)*mass1;
		
			//printf("masses=%e %e\nCM=(%e %e %e)\nmomentum=(%e %e %e)\n",mass1,mass2,cmpx,cmpy,cmpz,px,py,pz);
		
			// Find the collision orientation.
			dx=getx(p2)-getx(p1);
			dy=gety(p2)-gety(p1);
			dz=getz(p2)-getz(p1);
			dmag=sqrt(dx*dx+dy*dy+dz*dz);
			if(dmag==0.0) {
				printf("Error: Particles with zero separation!");
			} else if(dmag>1.1*(radius[p1.i]+radius[p2.i])) {
				printf("Error: We have a collision with separation of %e for particles of size %e and %e between %d and %d\n",dmag,radius[p1.i],radius[p2.i],p1.i,p2.i);
				return;
			} else if(dmag*1.0000001<radius[p1.i]+radius[p2.i]) { // overlapping so pull apart
//**				printf("Pulling apart %d and %d - %e + %e > %e\n",p1,p2,radius[p1.i],radius[p2.i],dmag);
				double move=(radius[p1.i]+radius[p2.i]-dmag)*0.5; // half the overlap
				dx=(getx(p2)-getx(p1))/dmag;
				dy=(gety(p2)-gety(p1))/dmag;
				dz=(getz(p2)-getz(p1))/dmag;
				cart[p1.i].p[0]-=dx*move;
				cart[p1.i].p[1]-=dy*move;
				cart[p1.i].p[2]-=dz*move;
				cart[p2.i].p[0]+=dx*move;
				cart[p2.i].p[1]+=dy*move;
				cart[p2.i].p[2]+=dz*move;
				dx=getx(p2)-getx(p1);
				dy=gety(p2)-gety(p1);
				dz=getz(p2)-getz(p1);
				dmag=sqrt(dx*dx+dy*dy+dz*dz);
			}
			dx=(getx(p2)-getx(p1))/dmag;
			dy=(gety(p2)-gety(p1))/dmag;
			dz=(getz(p2)-getz(p1))/dmag;
			//deltax[0]=fabs(getx(p2)-getx(p1));
			//deltax[1]=fabs(gety(p2)-gety(p1));
			//deltax[2]=fabs(getz(p2)-getz(p1));
			//printf("Sep %e %e %e\n",dx,dy,dz);
		
			// Check if this is a bad collision.
			gn=dx*px+dy*py+dz*pz;
			if(gn<0.0) {
//**				printf("Thrown out for wrong direction: %d %d gn=%e\n",p1,p2,gn/sqrt(px*px+py*py+pz*pz));
		//		printf("%e: %e %e  - %e %e %e\n",dmag,radius[p1.i],radius[p2.i],dx,dy,dz);
		//		printf("%e %e %e\n",px,py,pz);
				time[p1.i]=t;
				time[p2.i]=t;
				return;
			}

			// Calculate the new CM velocity.
			gnx=gn*dx;
			gny=gn*dy;
			gnz=gn*dz;
			gtx=px-gnx;
			gty=py-gny;
			gtz=pz-gnz;
			double velCmPerSec = fabs(gn/mass1+gn/mass2)*VelocityToCMperS();
		//	printf("gn=%e - (%e %e %e)\n",gn,gnx,gny,gnz);
			if(dmag>0.95*(radius[p1.i]+radius[p2.i])) {
				epsilon=CoefRest::EpsilonN(velCmPerSec);
			} else epsilon=1.0;
			gnx*=-epsilon;
			gny*=-epsilon;
			gnz*=-epsilon;
			px=gnx+gtx;
			py=gny+gty;
			pz=gnz+gtz;
		//	printf("epsilon=%e result P=(%e %e %e)\n",epsilon,px,py,pz);
		
			// Calculate the final velocity of the particles.
			cart[p1.i].p[3]=cmpx+px/mass1;
			cart[p1.i].p[4]=cmpy+py/mass1;
			cart[p1.i].p[5]=cmpz+pz/mass1;
			cart[p2.i].p[3]=cmpx-px/mass2;
			cart[p2.i].p[4]=cmpy-py/mass2;
			cart[p2.i].p[5]=cmpz-pz/mass2;
		//	printf("Final vels (%e %e %e)\n   (%e %e %e)\n\n",cart[p1.i].vx,cart[p1.i].vy,cart[p1.i].vz,cart[p2.i].vx,cart[p2.i].vy,cart[p2.i].vz);
			if((getvx(p2)-getvx(p1))*(getx(p2)-getx(p1))+(getvy(p2)-getvy(p1))*(gety(p2)-gety(p1))+(getvz(p2)-getvz(p1))*(getz(p2)-getz(p1))<0.0) {
				printf("After collision particles heading toward one another %d %d.\n",p1.i,p2.i);
				printf("Locations p1 - %e %e %e %e %e %e\n",getx(p1),gety(p1),getz(p1),getvx(p1),getvy(p1),getvz(p1));
				printf("Locations p2 - %e %e %e %e %e %e\n",getx(p2),gety(p2),getz(p2),getvx(p2),getvy(p2),getvz(p2));
			}
		
		//	printf("3. Average X=%1.16e\n",0.5*(gc[p1.i].X+gc[p2.i].X));
		}

		int getNumBodies() const { return numBodies; }
		void setNumBodies(int num) { 
			numBodies=num;
			cart.resize(num);
			radius.resize(num);
			time.resize(num);
		}
		int getNumReal() const { return numReal; }
		void setNumReal(int num) { numReal=num; }

		double get(ParticleIndex pi,int d) { return cart[pi.i].p[d]; }
		double getx(ParticleIndex pi) const { return cart[pi.i].p[0]; }
		double gety(ParticleIndex pi) const { return cart[pi.i].p[1]; }
		double getz(ParticleIndex pi) const { return cart[pi.i].p[2]; }
		double getvx(ParticleIndex pi) const { return cart[pi.i].p[3]; }
		double getvy(ParticleIndex pi) const { return cart[pi.i].p[4]; }
		double getvz(ParticleIndex pi) const { return cart[pi.i].p[5]; }

		void setx(ParticleIndex pi,double nx) { cart[pi.i].p[0]=nx; }
		void sety(ParticleIndex pi,double ny) { cart[pi.i].p[1]=ny; }
		void setz(ParticleIndex pi,double nz) { cart[pi.i].p[2]=nz; }
		void setvx(ParticleIndex pi,double nvx) { cart[pi.i].p[3]=nvx; }
		void setvy(ParticleIndex pi,double nvy) { cart[pi.i].p[4]=nvy; }
		void setvz(ParticleIndex pi,double nvz) { cart[pi.i].p[5]=nvz; }

		double getRadius(ParticleIndex pi) const { return radius[pi.i]; }
		void setRadius(ParticleIndex pi,double r) { radius[pi.i]=r; }
		double getTime(ParticleIndex pi) const { return time[pi.i]; }
		double getMass(ParticleIndex pi) const {
			return massFunc(pi, getRadius(pi));
		}

		double getTimeStep() const { return dt; }

		double getMaxParticleRadius() const { return maxRadius; }
		double getMinx() const { return bounds.getMinX(); }
		double getMaxx() const { return bounds.getMaxX(); }
		double getMiny() const { return bounds.getMinY(); }
		double getMaxy() const { return bounds.getMaxY(); }
		double getMinz() const { return bounds.getMinZ(); }
		double getMaxz() const { return bounds.getMaxZ(); }
		
		// Functions I added so this would work with Parallel Boundary
		void adjustAfterForce(ParticleIndex part) {}
		void setPhi(int i, double np){}
		void setZeta(int i, double nk){}
		void sete(int i, double ni){}
		void seti(int i, double ni){}
		void setY(int i, double ny){}
		void setX(int i, double nX){}
		double getX(int i) { return 0;}
		double getY(int i) {return 0;}
		double geti(int i) {return 0;}
		double gete(int i) {return 0;}
		double getZeta(int i) {return 0;}
		double getPhi(int i) {return 0;}
 
		std::vector<BasicCartCoords> &getCart() { return cart; }
		std::vector<double> &getRadius() { return radius; }

		void advanceParticleTo(ParticleIndex pi,double t) {
			AdvanceParticle(t-time[pi.i],pi);
		}

		double VelocityToCMperS() { return velTocms; }

	private:
		void AdvanceParticle(double t,ParticleIndex pi) {
			cart[pi.i].p[0]+=t*cart[pi.i].p[3];
			cart[pi.i].p[1]+=t*cart[pi.i].p[4];
#ifdef CONST_ACCEL
			cart[pi.i].p[5]+=t*CONST_ACCEL;
#endif
			cart[pi.i].p[2]+=t*cart[pi.i].p[5];
			time[pi.i]+=t;
		}

		int numBodies;
		int numReal;
		std::vector<BasicCartCoords> cart;
		std::vector<double> radius;
		std::vector<double> time;
		BoundaryCondition &bounds;
		OutputMethod &out;
		double dt;
		double velTocms;
		MassFunc massFunc;

		double maxRadius;
		
		long checkCnt;
};

#endif

