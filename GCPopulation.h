// GCPopulation.h
// This file describes a class that has a population of particles whose motion
// is described using guiding center variables.

#ifndef GC_POPULATION
#define GC_POPULATION

#include <vector>
#include <algorithm>
#include "Coordinates.h"
#include "CoefficientsOfRestitution.h"
#include "ParticleIndex.h"
#include "MassFunctions.h"

// If you want to have a file with the information on high velocity impacts
// you should define HIGH_VEL_OUTPUT and give it a numeric value of the
// number of collisions to output for each step.
#ifdef HIGH_VEL_OUTPUT
#include<queue>

struct CollVelData {
	double vel;
	double x1, y1, z1, vx1, vy1, vz1, rad1;
	double x2, y2, z2, vx2, vy2, vz2, rad2;
	bool operator<(const CollVelData &that) const {
		return vel > that.vel;
	}
};
#endif
		

template<class BoundaryCondition,class CollisionFinder,class OutputMethod,class GCType=GCCoords,class CoefRest=En_Bridges_Et_0p5, class MassFunc=StandardMass>
class GCPopulation {
	public:
		GCPopulation(BoundaryCondition &bc,OutputMethod &om,double timeStep,double r, MassFunc mf, double centralMass=5.68e26):
				bounds(bc),out(om),dt(timeStep),r0(r),massFunc(mf) {
			numBodies=0;
			numReal=0;
			maxRadius=0.0;
			
			r0*=1000;
			velTocms=100*sqrt(6.67e-11*centralMass/r0);
#ifdef DEBUG
			printf("%e, %e, %e\n", centralMass, r0, velTocms);
#endif
			
#ifdef HIGH_VEL_OUTPUT
			collVelFout = fopen("HighVelColls.bin", "wb");
#endif

			omp_init_lock(&lock);
		}

		~GCPopulation() {
#ifdef HIGH_VEL_OUTPUT
			fclose(collVelFout);
#endif
		}


		// For this functions the Distribution class needs to have two basic
		// pieces of functionality.  It needs to have a function that can
		// create a new particle as well as a function that can tell if another
		// particle needs to be generated.  This latter one allows simulations
		// that aim for a certain average optical depth, but it makes dealing
		// with the arrays more difficult.
		template<class Distribution>
		void randomDistribution(Distribution &pd) {
//			numBodies=0;
			double area=0.0;
			while(pd.moreParticles()) {
				numBodies++;
				gc.resize(numBodies);
				cart.resize(numBodies);
				radius.resize(numBodies);
				time.resize(numBodies);
#ifdef SPIN
				omega.resize(numBodies);
#endif
				if(pd.usesGC()) {
					pd.setNextParticle(gc[numBodies-1],radius[numBodies-1]);
					cart[numBodies-1].set(gc[numBodies-1]);
				} else {
#ifdef SPIN
					pd.setNextParticle(cart[numBodies-1],radius[numBodies-1],omega[numBodies-1]);
#else
					pd.setNextParticle(cart[numBodies-1],radius[numBodies-1]);
#endif
					gc[numBodies-1].set(cart[numBodies-1]);
				}
				area+=3.14159*radius[numBodies-1]*radius[numBodies-1];
				if(numBodies==1 || radius[numBodies-1]>maxRadius)
					maxRadius=radius[numBodies-1];
			}
			printf("Optical depth is %e\n",area/((getMaxx()-getMinx())*(getMaxy()-getMiny())));
			numReal=numBodies;
		}
		
		void addSingleParticleGC(double X,double Y,double e,double i,double phi,double zeta,double rad) {
			numBodies++;
			numReal++;
			gc.resize(numBodies);
			cart.resize(numBodies);
			radius.resize(numBodies);
			time.resize(numBodies);
			gc[numBodies-1].X=X;
			gc[numBodies-1].Y=Y;
			gc[numBodies-1].e=e;
			gc[numBodies-1].i=i;
			gc[numBodies-1].phi=phi;
			gc[numBodies-1].zeta=zeta;
			radius[numBodies-1]=rad;
			cart[numBodies-1].set(gc[numBodies-1]);
			if(numBodies==1 || radius[numBodies-1]>maxRadius)
				maxRadius=radius[numBodies-1];
#ifdef SPIN
			omega.resize(numBodies);
#endif
		}

		void addSingleParticleCart(double x,double y,double z,double vx,double vy,double vz,double rad) {
			numBodies++;
			numReal++;
			gc.resize(numBodies);
			cart.resize(numBodies);
			radius.resize(numBodies);
			time.resize(numBodies);
			cart[numBodies-1].x=x;
			cart[numBodies-1].y=y;
			cart[numBodies-1].z=z;
			cart[numBodies-1].vx=vx;
			cart[numBodies-1].vy=vy;
			cart[numBodies-1].vz=vz;
			radius[numBodies-1]=rad;
			gc[numBodies-1].set(cart[numBodies-1]);
			if(numBodies==1 || radius[numBodies-1]>maxRadius)
				maxRadius=radius[numBodies-1];
#ifdef SPIN
			omega.resize(numBodies);
#endif
		}

#ifdef SPIN
		void addSingleParticleCartSpin(double x,double y,double z,double vx,double vy,double vz,double rad,double wx,double wy,double wz) {
			numBodies++;
			numReal++;
			gc.resize(numBodies);
			cart.resize(numBodies);
			radius.resize(numBodies);
			time.resize(numBodies);
			cart[numBodies-1].x=x;
			cart[numBodies-1].y=y;
			cart[numBodies-1].z=z;
			cart[numBodies-1].vx=vx;
			cart[numBodies-1].vy=vy;
			cart[numBodies-1].vz=vz;
			radius[numBodies-1]=rad;
			gc[numBodies-1].set(cart[numBodies-1]);
			if(numBodies==1 || radius[numBodies-1]>maxRadius)
				maxRadius=radius[numBodies-1];
			omega.resize(numBodies);
			omega[numBodies-1].x=wx;
			omega[numBodies-1].y=wy;
			omega[numBodies-1].z=wz;
		}
#endif

		/**
		 * This method will add in a single particle then search the
		 * rest of the population and remove all the particles it
		 * overlapped with.
		 */
		void addSingleParticleGCandClear(double X,double Y,double e,double i,double phi,double zeta,double rad) {
			GCType newGC(X,Y,e,phi,i,zeta);
			CartCoords newCart(newGC);
			for(int j=0; j<numBodies; ++j) {
				if(newCart.dist(cart[j])<rad+radius[j]) {
					removeParticle(j);
					--j;
				}
			}
			addSingleParticleGC(X,Y,e,i,phi,zeta,rad);
		}
		
		void removeSphereOfParticlesCart(double x,double y,double z,double rad) {
			CartCoords newCart(x,y,z,0,0,0);
			for(int j=0; j<numBodies; ++j) {
				if(newCart.dist(cart[j])<rad+radius[j]) {
					removeParticle(j);
					--j;
				}
			}
		}

		void addLatticeParticleandClear(double x,double y,double z,double R,double r,double vx=0.0,double vy=0.0,double vz=0.0) {
			for(int j=0; j<numBodies; ++j) {
				double dx=x-cart[j].x;
				double dy=y-cart[j].y;
				double dz=z-cart[j].z;
				if(sqrt(dx*dx+dy*dy+dz*dz)<R+radius[j]) {
					removeParticle(j);
					--j;
				}
			}
			//printf("Lattice at %e %e %e with %e of %e\n",x,y,z,R,r);
			int num=addLatticeSphere(x,y,z,R,r,vx,vy,vz);
			printf("Added Lattice of %d particles.\n",num);
		}

		void addLatticeParticleandClear(double x,double y,double z,double Rx,double Ry,double Rz,double r,double vx=0.0,double vy=0.0,double vz=0.0) {
			for(int j=0; j<numBodies; ++j) {
				double dx=x-cart[j].x;
				double dy=y-cart[j].y;
				double dz=z-cart[j].z;
				if(sqrt(dx*dx+dy*dy*Rx/Ry+dz*dz*Rx/Rz)<Rx+radius[j]) {
					removeParticle(j);
					--j;
				}
			}
			//printf("Lattice at %e %e %e with %e of %e\n",x,y,z,R,r);
			int num=addLatticeSphere(x,y,z,Rx,Ry,Rz,r,vx,vy,vz);
			printf("Added Lattice of %d particles.\n",num);
		}

		int addLatticeSphere(double x,double y,double z,double R,double r,double vx=0.0,double vy=0.0,double vz=0.0) {
			double xoff=1.0001*r;
			double yoff=1.0001*r*sqrt(3)/3;
			double zoff=1.0001*2*r*sqrt(6)/3;
			double fx=x,fy=y,fz=z;
			int n1=100,n2=100;
			int ret=0;

			while(n1>0 && n2>0) {
				n1=addLatticePlane(fx,fy,fz,x,y,z,R,R,R,r,vx,vy,vz);
				n2=addLatticePlane(fx+xoff,fy+yoff,fz+zoff,x,y,z,R,R,R,r,vx,vy,vz);
				fz+=2*zoff;
				ret+=n1+n2;
			}
			fz=z-2*zoff;
			n1=n2=100;
			while(n1>0 && n2>0) {
				n1=addLatticePlane(fx,fy,fz,x,y,z,R,R,R,r,vx,vy,vz);
				n2=addLatticePlane(fx+xoff,fy+yoff,fz+zoff,x,y,z,R,R,R,r,vx,vy,vz);
				fz-=2*zoff;
				ret+=n1+n2;
			}
			return ret;
		}

		int addLatticeEllipsoid(double x,double y,double z,double Rx,double Ry,double Rz,double r,double vx=0.0,double vy=0.0,double vz=0.0) {
			double xoff=1.0001*r;
			double yoff=1.0001*r*sqrt(3)/3;
			double zoff=1.0001*2*r*sqrt(6)/3;
			double fx=x,fy=y,fz=z;
			int n1=100,n2=100;
			int ret=0;

			while(n1>0 && n2>0) {
				n1=addLatticePlane(fx,fy,fz,x,y,z,Rx,Ry,Rz,r,vx,vy,vz);
				n2=addLatticePlane(fx+xoff,fy+yoff,fz+zoff,x,y,z,Rx,Ry,Rz,r,vx,vy,vz);
				fz+=2*zoff;
				ret+=n1+n2;
			}
			fz=z-2*zoff;
			n1=n2=100;
			while(n1>0 && n2>0) {
				n1=addLatticePlane(fx,fy,fz,x,y,z,Rx,Ry,Rz,r,vx,vy,vz);
				n2=addLatticePlane(fx+xoff,fy+yoff,fz+zoff,x,y,z,Rx,Ry,Rz,r,vx,vy,vz);
				fz-=2*zoff;
				ret+=n1+n2;
			}
			return ret;
		}

		int addLatticePlane(double fx,double fy,double fz,double x,double y,double z,double Rx,double Ry,double Rz,double r,double vx,double vy,double vz) {
			int ret=0;
			double yoff=1.0001*r*sqrt(3);
			double xoff=1.0001*r;
			int n1=100,n2=100;
			double fy_i=fy;

			while(n1>0 && n2>0) {
				n1=addLatticeLine(fx,fy,fz,x,y,z,Rx,Ry,Rz,r,vx,vy,vz);
				n2=addLatticeLine(fx+xoff,fy+yoff,fz,x,y,z,Rx,Ry,Rz,r,vx,vy,vz);
				ret+=n1+n2;
				fy+=2*yoff;
			}
			fy=fy_i-2*yoff;
			n1=n2=100;
			while(n1>0 && n2>0) {
				n1=addLatticeLine(fx,fy,fz,x,y,z,Rx,Ry,Rz,r,vx,vy,vz);
				n2=addLatticeLine(fx+xoff,fy+yoff,fz,x,y,z,Rx,Ry,Rz,r,vx,vy,vz);
				ret+=n1+n2;
				fy-=2*yoff;
			}
			return ret;
		}

		int addLatticeLine(double fx,double fy,double fz,double x,double y,double z,double Rx,double Ry,double Rz,double r,double vx,double vy,double vz) {
			int ret=0;
			double fx_i=fx;
			double dx=fx-x;
			double dy=fy-y;
			double dz=fz-z;
			//printf("Lattice Start line %e %e %e\n",fx,fy,fz);
			while(sqrt(dx*dx+dy*dy*Rx/Ry*Rx/Ry+dz*dz*Rx/Rz*Rx/Rz)<Rx-r) {
				addSingleParticleCart(fx,fy,fz,vx,vy,vz,r);
				//printf("Lattice Adding at %e %e %e\n",fx,fy,fz);
				ret++;
				fx+=1.0001*2*r;
				dx=fx-x;
			}
			fx=fx_i-1.0001*2*r;
			dx=fx-x;
			while(sqrt(dx*dx+dy*dy*Rx/Ry*Rx/Ry+dz*dz*Rx/Rz*Rx/Rz)<Rx-r) {
				addSingleParticleCart(fx,fy,fz,vx,vy,vz,r);
				//printf("Lattice Adding at %e %e %e\n",fx,fy,fz);
				ret++;
				fx-=1.0001*2*r;
				dx=fx-x;
			}
			return ret;
		}
		
		/**
		 * This function will remove a single particle and take the
		 * particle from the end and move it into that place.  You
		 * can't do this safely in the middle of a timestep as it
		 * changes the indexes of particles.
		 */
		void removeParticle(int index) {
			omp_set_lock(&lock);
			numBodies--;
			numReal--;
			gc[index]=gc[numReal];
			cart[index]=cart[numReal];
			radius[index]=radius[numReal];
			time[index]=time[numReal];
#ifdef SPIN
			omega[index]=omega[numReal];
#endif
			if(numBodies>numReal) {
				gc[numReal]=gc[numBodies];
				cart[numReal]=cart[numBodies];
				radius[numReal]=radius[numBodies];
				time[numReal]=time[numBodies];
#ifdef SPIN
				omega[numReal]=omega[numBodies];
#endif
			}
			gc.resize(gc.size()-1);
			cart.resize(cart.size()-1);
			radius.resize(radius.size()-1);
			time.resize(time.size()-1);
			omp_unset_lock(&lock);
		}
		
		// Right now I'm not using an endCart array.  This means more
		// coordinate conversions for doing collisions to full accuracy.
		void endStep(void) {
			// Advance all particles to end of step.
			#pragma omp parallel for schedule(static)
			for(int i=0; i<numBodies; ++i) {
				double t=dt-time[i];
				gc[i].advance(t);
				while(gc[i].phi>6.283185307) gc[i].phi-=6.283185307;
				while(gc[i].zeta>6.283185307) gc[i].zeta-=6.283185307;
				cart[i].set(gc[i]);
				time[i]=0.0;
			}

			bounds.apply(*this);
			out.output(*this);

#ifdef HIGH_VEL_OUTPUT
			int s = collPQ.size();
			fwrite(&s, sizeof(int), 1, collVelFout);
			while(!collPQ.empty()) {
				fwrite(&(collPQ.top()), sizeof(CollVelData), 1, collVelFout);
				collPQ.pop();
			}
			fflush(collVelFout);
#endif
		}

		// This adjusts important parameters for the given particle based
		// on the new cartesian coordinates.
		void adjustAfterForce(ParticleIndex part) {
			gc[part.i].set(cart[part.i]);
		}

		// This adjusts important parameters for the given particle based
		// on the new guiding center coordinates.
		void setCartAfterForce(ParticleIndex part) {
			cart[part.i].set(gc[part.i]);
		}

		double collisionTime(ParticleIndex p1,ParticleIndex p2) {
			if(time[p1.i]==0 && time[p2.i]==0)
				return CollisionFinder::collisionTimeInitial(cart[p1.i],cart[p2.i],gc[p1.i],gc[p2.i],radius[p1.i],radius[p2.i],dt);
			else
				return CollisionFinder::collisionTime(cart[p1.i],cart[p2.i],gc[p1.i],gc[p2.i],radius[p1.i],radius[p2.i],time[p1.i],time[p2.i],dt);
		}

		void processCollision(ParticleIndex p1,ParticleIndex p2,double t) {
			double cmpx,cmpy,cmpz,dx,dy,dz,gn,epsilon;
			double gnx,gny,gnz,gtx,gty,gtz;
			double px,py,pz;
			double dmag;
			double mass1=getMass(p1);
			double mass2=getMass(p2);
			double tmass=mass1+mass2;
			//double deltax[3];
		
		//	printf("1. Average X=%1.16e\n",0.5*(gc[p1].X+gc[p2].X));
//			printf("Process %e/%e, %d, %d\n",t,dt,p1,p2);
		
			// Advance the particles to the time of the collision.
			if((t<time[p1.i]) || (t<time[p2.i])) {
				printf("%d [%d,%d] at %e < %e or %e: Bad time in collision!\n",omp_get_thread_num(),p1.i,p2.i,t,time[p1.i],time[p2.i]);
				return;
			}
		
			AdvanceParticle(t-time[p1.i],p1);
			AdvanceParticle(t-time[p2.i],p2);
		
		//	printf("Collision between %d and %d at %e\n",p1,p2,t);
		//	printf("Initial vels (%e %e %e)\n   (%e %e %e)\n",cart[p1].vx,cart[p1].vy,cart[p1].vz,cart[p2].vx,cart[p2].vy,cart[p2].vz);
			// Move to the center of mass frame.
			cmpx=(mass2*cart[p2.i].vx+mass1*cart[p1.i].vx)/tmass;
			cmpy=(mass2*cart[p2.i].vy+mass1*cart[p1.i].vy)/tmass;
			cmpz=(mass2*cart[p2.i].vz+mass1*cart[p1.i].vz)/tmass;
			px=(cart[p1.i].vx-cmpx)*mass1;
			py=(cart[p1.i].vy-cmpy)*mass1;
			pz=(cart[p1.i].vz-cmpz)*mass1;

#ifdef HIGH_VEL_OUTPUT
			{  // Create this scope so that variables don't conflict
				double vx = cart[p2.i].vx-cart[p1.i].vx;
				double vy = cart[p2.i].vy-cart[p1.i].vy;
				double vz = cart[p2.i].vz-cart[p1.i].vz;
				double vel = sqrt(vx*vx + vy*vy + vz*vz);
				omp_set_lock(&lock);
				if(collPQ.size() < HIGH_VEL_OUTPUT || vel > collPQ.top().vel) {
					collPQ.push(
						CollVelData {	vel,
							cart[p1.i].x, cart[p1.i].y, cart[p1.i].z,
							cart[p1.i].vx, cart[p1.i].vy, cart[p1.i].vz,
							radius[p1.i],
							cart[p2.i].x, cart[p2.i].y, cart[p2.i].z,
							cart[p2.i].vx, cart[p2.i].vy, cart[p2.i].vz,
							radius[p2.i]
						});
					if(collPQ.size() > HIGH_VEL_OUTPUT) collPQ.pop();
				}
				omp_unset_lock(&lock);
			}
#endif
		
			//printf("masses=%e %e\nCM=(%e %e %e)\nmomentum=(%e %e %e)\n",mass1,mass2,cmpx,cmpy,cmpz,px,py,pz);
		
			// Find the collision orientation.
			dx=cart[p1.i].x-cart[p2.i].x;
			dy=cart[p1.i].y-cart[p2.i].y;
			dz=cart[p1.i].z-cart[p2.i].z;
			dmag=sqrt(dx*dx+dy*dy+dz*dz);
			if(dmag > 1.1*(radius[p1.i]+radius[p2.i])) {
				printf("%d [%d, %d] at %e: We have a collision with separation of %e for particles of size %e and %e\n",omp_get_thread_num(),p1.i,p2.i,t,dmag,radius[p1.i],radius[p2.i]);
				return;
			} else if(dmag < radius[p1.i]+radius[p2.i]) { // overlapping so pull apart
#ifdef DEBUG
				if (dmag < (radius[p1.i]+radius[p2.i])*0.9) {
					printf("ERROR!! Serious overlap. %d %d %f\n", p1.i, p2.i, dmag / (radius[p1.i]+radius[p2.i]));
				}
#endif
				//printf("Overlap: pushing apart %e<%e+%e\n",dmag,radius[p1],radius[p2]);
				double move=radius[p1.i]+radius[p2.i]-dmag; // half the overlap
				double move1=mass2/tmass;
				double move2=1.0-move1;
				move1*=move;
				move2*=move;
				dx=(cart[p2.i].x-cart[p1.i].x)/dmag;
				dy=(cart[p2.i].y-cart[p1.i].y)/dmag;
				dz=(cart[p2.i].z-cart[p1.i].z)/dmag;
				cart[p1.i].x-=dx*move1;
				cart[p1.i].y-=dy*move1;
				cart[p1.i].z-=dz*move1;
				cart[p2.i].x+=dx*move2;
				cart[p2.i].y+=dy*move2;
				cart[p2.i].z+=dz*move2;
				// Set the carts here in case gn < 0 and we jump out early.
				gc[p1.i].set(cart[p1.i]);
				gc[p2.i].set(cart[p2.i]);
				dx=cart[p1.i].x-cart[p2.i].x;
				dy=cart[p1.i].y-cart[p2.i].y;
				dz=cart[p1.i].z-cart[p2.i].z;
				dmag=sqrt(dx*dx+dy*dy+dz*dz);
#ifdef DEBUG
				if (dmag < (radius[p1.i]+radius[p2.i]*0.999)) {
					printf("ERROR!! Overlap after push apart. %d %d\n", p1.i, p2.i);
				}
#endif
			}
			dx=(cart[p2.i].x-cart[p1.i].x)/dmag;
			dy=(cart[p2.i].y-cart[p1.i].y)/dmag;
			dz=(cart[p2.i].z-cart[p1.i].z)/dmag;
			//deltax[0]=fabs(cart[p2].x-cart[p1].x);
			//deltax[1]=fabs(cart[p2].y-cart[p1].y);
			//deltax[2]=fabs(cart[p2].z-cart[p1].z);
			//printf("Sep %e %e %e\n",dx,dy,dz);
		
			// Check if this is a bad collision.
			gn=dx*px+dy*py+dz*pz;
			if(gn<0.0) {
				double magp=sqrt(px*px+py*py+pz*pz);
				if(gn/magp<-1e-3)
					printf("%d [%d,%d] at %e with %e %e: Thrown out for wrong direction: gn=%e\n",omp_get_thread_num(),p1.i,p2.i,t,time[p1.i],time[p2.i],gn/magp);
		//		printf("%e: %e %e  - %e %e %e\n",dmag,radius[p1],radius[p2],dx,dy,dz);
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
			epsilon=CoefRest::EpsilonN(velCmPerSec);
#ifdef DEBUG
			printf("vel = %e, epsilon = %e\n", velCmPerSec, epsilon);
			printf("convert = %e\n", VelocityToCMperS());
#endif
		//		if(mass1>mass2) epsilon*=mass2/mass1;
		//		else if(mass2>mass1) epsilon*=mass1/mass2;
//			if(radius[p1]>2e-7 || radius[p2]>2e-7) printf("epsilon=%e gn=%e result P=(%e %e %e)\n",epsilon,gn,px,py,pz);
#ifdef SPIN
			double vx = cart[p2.i].vx-cart[p1.i].vx;
			double vy = cart[p2.i].vy-cart[p1.i].vy;
			double vz = cart[p2.i].vz-cart[p1.i].vz;
			double R1x = radius[p1.i]*dx;
			double R1y = radius[p1.i]*dy;
			double R1z = radius[p1.i]*dz;
			double sigma1x = omega[p1.i].y*R1z - omega[p1.i].z*R1y;
			double sigma1y = omega[p1.i].z*R1x - omega[p1.i].x*R1z;
			double sigma1z = omega[p1.i].x*R1y - omega[p1.i].y*R1x;
			double R2x = -radius[p2.i]*dx;
			double R2y = -radius[p2.i]*dy;
			double R2z = -radius[p2.i]*dz;
			double sigma2x = omega[p2.i].y*R2z - omega[p2.i].z*R2y;
			double sigma2y = omega[p2.i].z*R2x - omega[p2.i].x*R2z;
			double sigma2z = omega[p2.i].x*R2y - omega[p2.i].y*R2x;
			double sigmax = sigma2x-sigma1x;
			double sigmay = sigma2y-sigma1y;
			double sigmaz = sigma2z-sigma1z;
			double ux = vx+sigmax;
			double uy = vy+sigmay;
			double uz = vz+sigmaz;
			double udotn = ux*dx+uy*dy+uz*dz;
			double unx = dx*udotn;
			double uny = dy*udotn;
			double unz = dz*udotn;
			double utx = ux-unx;
			double uty = uy-uny;
			double utz = uz-unz;

#ifdef DEBUG
			printf("omega1 = %e, %e, %e\n",omega[p1.i].x,omega[p1.i].y,omega[p1.i].z);
			printf("omega2 = %e, %e, %e\n",omega[p2.i].x,omega[p2.i].y,omega[p2.i].z);
			printf("v = %e, %e, %e\n",vx,vy,vz);
			printf("R1 = %e, %e, %e\n",R1x,R1y,R1z);
			printf("sigma1 = %e, %e, %e\n",sigma1x,sigma1y,sigma1z);
			printf("R2 = %e, %e, %e\n",R2x,R2y,R2z);
			printf("sigma2 = %e, %e, %e\n",sigma2x,sigma2y,sigma2z);
			printf("sigma = %e, %e, %e\n",sigmax,sigmay,sigmaz);
			printf("u = %e, %e, %e\n",ux,uy,uz);
			printf("un = %e, %e, %e\n",unx,uny,unz);
			printf("ut = %e, %e, %e\n",utx,uty,utz);
#endif

			double epsilon_t = CoefRest::EpsilonT(sqrt(utx * utx +
				uty * uty + utz * utz)*VelocityToCMperS());
			double beta = 2.0/7.0;

			cart[p1.i].vx += mass2/tmass*((1+epsilon)*unx+beta*(1-epsilon_t)*utx);
			cart[p1.i].vy += mass2/tmass*((1+epsilon)*uny+beta*(1-epsilon_t)*uty);
			cart[p1.i].vz += mass2/tmass*((1+epsilon)*unz+beta*(1-epsilon_t)*utz);
			cart[p2.i].vx -= mass1/tmass*((1+epsilon)*unx+beta*(1-epsilon_t)*utx);
			cart[p2.i].vy -= mass1/tmass*((1+epsilon)*uny+beta*(1-epsilon_t)*uty);
			cart[p2.i].vz -= mass1/tmass*((1+epsilon)*unz+beta*(1-epsilon_t)*utz);

			double mu = mass1*mass2/tmass;
			double I1 = 0.4*mass1*radius[p1.i]*radius[p1.i];
			double I2 = 0.4*mass2*radius[p2.i]*radius[p2.i];
			omega[p1.i].x += beta*mu/I1*(1-epsilon_t)*(R1y*uz-R1z*uy);
			omega[p1.i].y += beta*mu/I1*(1-epsilon_t)*(R1z*ux-R1x*uz);
			omega[p1.i].z += beta*mu/I1*(1-epsilon_t)*(R1x*uy-R1y*ux);
			omega[p2.i].x -= beta*mu/I2*(1-epsilon_t)*(R2y*uz-R2z*uy);
			omega[p2.i].y -= beta*mu/I2*(1-epsilon_t)*(R2z*ux-R2x*uz);
			omega[p2.i].z -= beta*mu/I2*(1-epsilon_t)*(R2x*uy-R2y*ux);
//			printf("omega1 = %e, %e, %e\n",omega[p1].x,omega[p1].y,omega[p1].z);
//			printf("omega2 = %e, %e, %e\n",omega[p2].x,omega[p2].y,omega[p2].z);
#else
			gnx*=-epsilon;
			gny*=-epsilon;
			gnz*=-epsilon;
			px=gnx+gtx;
			py=gny+gty;
			pz=gnz+gtz;
			// Calculate the final velocity of the particles.
			cart[p1.i].vx=cmpx+px/mass1;
			cart[p1.i].vy=cmpy+py/mass1;
			cart[p1.i].vz=cmpz+pz/mass1;
			cart[p2.i].vx=cmpx-px/mass2;
			cart[p2.i].vy=cmpy-py/mass2;
			cart[p2.i].vz=cmpz-pz/mass2;
#endif

#ifdef DEBUG
		//	printf("Final vels (%e %e %e)\n   (%e %e %e)\n\n",cart[p1].vx,cart[p1].vy,cart[p1].vz,cart[p2].vx,cart[p2].vy,cart[p2].vz);
			if((cart[p2.i].vx-cart[p1.i].vx)*(cart[p2.i].x-cart[p1.i].x)+(cart[p2.i].vy-cart[p1.i].vy)*(cart[p2.i].y-cart[p1.i].y)+(cart[p2.i].vz-cart[p1.i].vz)*(cart[p2.i].z-cart[p1.i].z)<0.0) {
				printf("After collision particles heading toward one another %d %d.\n",p1.i,p2.i);
				printf("Locations p1 - %e %e %e %e %e %e\n",cart[p1.i].x,cart[p1.i].y,cart[p1.i].z,cart[p1.i].vx,cart[p1.i].vy,cart[p1.i].vz);
				printf("Locations p2 - %e %e %e %e %e %e\n",cart[p2.i].x,cart[p2.i].y,cart[p2.i].z,cart[p2.i].vx,cart[p2.i].vy,cart[p2.i].vz);
				if(p1.i>=5000 || p2.i>=5000) exit(1);
			}
			if(sqrt((cart[p1.i].x-cart[p2.i].x)*(cart[p1.i].x-cart[p2.i].x) + 
				 (cart[p1.i].y-cart[p2.i].y)*(cart[p1.i].y-cart[p2.i].y) +
				 (cart[p1.i].z-cart[p2.i].z)*(cart[p1.i].z-cart[p2.i].z)) < 0.99*(radius[p1.i]+radius[p2.i])) {
				printf("Overlap at end of collision %d %d\n", p1.i, p2.i);
			}
#endif		
			// Set the new end positions
			gc[p1.i].set(cart[p1.i]);
			gc[p2.i].set(cart[p2.i]);
		//	printf("3. Average X=%1.16e\n",0.5*(gc[p1].X+gc[p2].X));
		}

		int getNumBodies() const { return numBodies; }
		void setNumBodies(int num) { 
			numBodies=num;
			gc.resize(num);
			cart.resize(num);
			radius.resize(num);
			time.resize(num);
#ifdef SPIN
			omega.resize(num);
#endif
		}

		int getNumReal() const { return numReal; }
		void setNumReal(int num) { numReal=num; }

		double getx(ParticleIndex pi) const { return cart[pi.i].x; }
		double gety(ParticleIndex pi) const { return cart[pi.i].y; }
		double getz(ParticleIndex pi) const { return cart[pi.i].z; }
		double getvx(ParticleIndex pi) const { return cart[pi.i].vx; }
		double getvy(ParticleIndex pi) const { return cart[pi.i].vy; }
		double getvz(ParticleIndex pi) const { return cart[pi.i].vz; }
		double get(ParticleIndex pi,int dim) {
//			return cart[i].get(dim);
			switch(dim) {
				case 0: return cart[pi.i].x;
				case 1: return cart[pi.i].y;
				case 2: return cart[pi.i].z;
				case 3: return cart[pi.i].vx;
				case 4: return cart[pi.i].vy;
				case 5: return cart[pi.i].vz;
//				case 6: return cart[i].fill1;
//				case 7: return cart[i].fill1;
			}
			printf("Bad dimension!!\n");
			return 0.0;
		}

		void setx(ParticleIndex pi,double nx) { cart[pi.i].x=nx; }
		void sety(ParticleIndex pi,double ny) { cart[pi.i].y=ny; }
		void setz(ParticleIndex pi,double nz) { cart[pi.i].z=nz; }
		void setvx(ParticleIndex pi,double nvx) { cart[pi.i].vx=nvx; }
		void setvy(ParticleIndex pi,double nvy) { cart[pi.i].vy=nvy; }
		void setvz(ParticleIndex pi,double nvz) { cart[pi.i].vz=nvz; }

		double getX(ParticleIndex pi) const { return gc[pi.i].X; }
		double getY(ParticleIndex pi) const { return gc[pi.i].Y; }
		double gete(ParticleIndex pi) const { return gc[pi.i].e; }
		double geti(ParticleIndex pi) const { return gc[pi.i].i; }
		double getPhi(ParticleIndex pi) const { return gc[pi.i].phi; }
		double getZeta(ParticleIndex pi) const { return gc[pi.i].zeta; }

		void setX(ParticleIndex pi,double nX) { gc[pi.i].X=nX; }
		void setY(ParticleIndex pi,double nY) { gc[pi.i].Y=nY; }
		void sete(ParticleIndex pi,double ne) { gc[pi.i].e=ne; }
		void seti(ParticleIndex pi,double ni) { gc[pi.i].i=ni; }
		void setPhi(ParticleIndex pi,double np) { gc[pi.i].phi=np; }
		void setZeta(ParticleIndex pi,double nk) { gc[pi.i].zeta=nk; }

		double getRadius(ParticleIndex pi) const {
			return radius[pi.i];
		}
		void setRadius(ParticleIndex pi,double r) { radius[pi.i]=r; }
		double getTime(ParticleIndex pi) const { return time[pi.i]; }
		double getMass(ParticleIndex pi) const {
			return massFunc(pi, getRadius(pi));
		}
#ifdef SPIN
		double getwx(ParticleIndex pi) { return omega[pi.i].x; }
		double getwy(ParticleIndex pi) { return omega[pi.i].y; }
		double getwz(ParticleIndex pi) { return omega[pi.i].z; }
		void setwx(ParticleIndex pi,double wx) { omega[pi.i].x=wx; }
		void setwy(ParticleIndex pi,double wy) { omega[pi.i].y=wy; }
		void setwz(ParticleIndex pi,double wz) { omega[pi.i].z=wz; }
#endif

		double getTimeStep() const { return dt; }

		double getMaxParticleRadius() const { return maxRadius; }
		double getMinx() const { return bounds.getMinX(); }
		double getMaxx() const { return bounds.getMaxX(); }
		double getMiny() const { return bounds.getMinY(); }
		double getMaxy() const { return bounds.getMaxY(); }

		std::vector<GCType> &getGC() { return gc; }
		std::vector<CartCoords> &getCart() { return cart; }
		std::vector<double> &getRadius() { return radius; }
#ifdef SPIN
		std::vector<SpinVector> &getomega() { return omega; }
#endif

		double VelocityToCMperS() { return velTocms; }

#ifdef REORDER
		void swapParticles(ParticleIndex pi,ParticleIndex pj) {
			std::swap(gc[pi.i],gc[pj.i]);
			std::swap(cart[pi.i],cart[pj.i]);
			std::swap(radius[pi.i],radius[pj.i]);
			std::swap(time[pi.i],time[pj.i]);
#ifdef SPIN
			std::swap(omega[pi.i],omega[pj.i]);
#endif
		}
#endif

	private:
		void AdvanceParticle(double t,ParticleIndex pi) {
			gc[pi.i].advance(t);
      while(gc[pi.i].phi>6.283185307) gc[pi.i].phi-=6.283185307;
      while(gc[pi.i].zeta>6.283185307) gc[pi.i].zeta-=6.283185307;
			cart[pi.i].set(gc[pi.i]);
			time[pi.i]+=t;
		}

		int numBodies;
		int numReal;
		std::vector<GCType> gc;
		std::vector<CartCoords> cart;
#ifdef SPIN
		std::vector<SpinVector> omega;
#endif
		std::vector<double> radius;
		std::vector<double> time;
		BoundaryCondition &bounds;
		OutputMethod &out;
#ifdef HIGH_VEL_OUTPUT
		std::priority_queue<CollVelData> collPQ;
		FILE *collVelFout;
#endif

		double dt;
		double maxRadius;
		double r0;
		double velTocms;
		MassFunc massFunc;

		omp_lock_t lock;
};

#endif

