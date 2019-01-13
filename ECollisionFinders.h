// CollisionFinders.h
// This file contains various classes that can be used to find collisions
// between particles.

#include <stdio.h>
#include <math.h>
#include "GCPopulation.h"

#ifndef ECOLLISION_FINDERS
#define ECOLLISION_FINDERS

class EFullNewtonFinder {
	public:
		static double collisionTime(const CartCoords &cart1,const CartCoords &cart2,const EGCCoords &gc1,const EGCCoords &gc2,double radius1,double radius2,double time1,double time2,double dt) {
			double max_time=(time1>time2)?(time1):(time2);
			double dist_sqr,dist_sqr_prime;
			double rad_sqr=(radius1+radius2)*(radius1+radius2);
			GCDistancePrime(gc1,gc2,max_time-time1,max_time-time2,&dist_sqr,&dist_sqr_prime);
			if(dist_sqr<rad_sqr) {
				if(dist_sqr_prime<0.0) {
					return(max_time+dt*0.0001);
				} else return(-2.0);
			}
		
			double guess=max_time;
//			GCDistancePrime(gc1,gc2,guess-time1,guess-time2,&dist_sqr,&dist_sqr_prime);
			if(dist_sqr_prime>0.0) return(-3.0);
			for(int i=0; i<50; i++) {
				if(fabs(dist_sqr-rad_sqr)<1e-8*rad_sqr) {
					if((dist_sqr_prime>0.0) && (guess>0.0)) {
						while(dist_sqr_prime>0.0) {
							guess=0.5*(guess+max_time);
							GCDistancePrime(gc1,gc2,guess-time1,guess-time2,&dist_sqr,&dist_sqr_prime);
						}
					} else {
		//				printf("Returning a guess of %e %e %e for %d %d\n",guess,dist_sqr,rad_sqr,p1,p2);
						return(guess);
					}
				}
				double step=-0.5*(dist_sqr-rad_sqr)/dist_sqr_prime;
				double guess2=guess+step;
				double dist_sqr2=GCDistance(gc1,gc2,guess2-time1,guess2-time2);
				int cnt=0;
				while((fabs(dist_sqr2-rad_sqr)>fabs(dist_sqr-rad_sqr)) && (cnt<15)) {
		//			printf("%d %e %e %e %e\n",i,step,guess,dist_sqr2-rad_sqr,dist_sqr-rad_sqr);
					step*=0.5;
					guess2=guess+step;
					dist_sqr2=GCDistance(gc1,gc2,guess2-time1,guess2-time2);
					cnt++;
				}
				if(cnt>=15) {
		//			printf("Cut in half too many times.\n");
					return(-7.0);
				}
				guess=guess2;
				GCDistancePrime(gc1,gc2,guess-time1,guess-time2,&dist_sqr,&dist_sqr_prime);
				if(guess<-5.0*dt) {
		//			printf("guess got too small %e\n",guess);
					return(-5.0);
				}
				if(guess>5.0*dt) {
		//			printf("guess got too large %e\n",guess);
					return(-6.0);
				}
			}
			printf("Didn't find a collision or jump out in 50 steps.  %e\n",guess);
			printf("Location 1 - %e %e %e %e %e %e\n",cart1.x,cart1.y,cart1.z,cart1.vx,cart1.vy,cart1.vz);
			printf("Location 2 - %e %e %e %e %e %e\n",cart2.x,cart2.y,cart2.z,cart2.vx,cart2.vy,cart2.vz);
			return(-1.0);
		}
	private:
		static double GCDistance(const EGCCoords &gc1,const EGCCoords &gc2,double t1,double t2) {
			double term1=1.0-1.5*gc1.X;
			double term2=1.0-1.5*gc2.X;
			double dx=gc1.X-gc2.X-gc1.e*cos(gc1.phi+t1*term1)+gc2.e*cos(gc2.phi+t2*term2);
			double dy=(gc1.Y-1.5*gc1.X*term1*t1)-(gc2.Y-1.5*gc2.X*term2*t2)+2.0*gc1.e*sin(gc1.phi+t1*term1)-2.0*gc2.e*sin(gc2.phi+t2*term2);
			double dz=gc1.i*cos(gc1.zeta+t1*term1)-gc2.i*cos(gc2.zeta+t2*term2);
			return(dx*dx+dy*dy+dz*dz);
		}

		static void GCDistancePrime(const EGCCoords &gc1,const EGCCoords &gc2,double t1,double t2,double *d,double *dp) {
			double term1=1.0-1.5*gc1.X;
			double term2=1.0-1.5*gc2.X;
			double cosPhi1=cos(gc1.phi+t1*term1);
			double cosPhi2=cos(gc2.phi+t2*term2);
			double sinPhi1=sin(gc1.phi+t1*term1);
			double sinPhi2=sin(gc2.phi+t2*term2);
			double dx=gc1.X-gc2.X-gc1.e*cosPhi1+gc2.e*cosPhi2;
			double dy=(gc1.Y-1.5*gc1.X*term1*t1)-(gc2.Y-1.5*gc2.X*term2*t2)+2.0*gc1.e*sinPhi1-2.0*gc2.e*sinPhi2;
			double dz=gc1.i*cos(gc1.zeta+t1*term1)-gc2.i*cos(gc2.zeta+t2*term2);
			*d=dx*dx+dy*dy+dz*dz;
			double dxp=term1*gc1.e*sinPhi1-term2*gc2.e*sinPhi2;
			double dyp=term1*(-1.5*gc1.X+2.0*gc1.e*cosPhi1)+term2*(1.5*gc2.X-2.0*gc2.e*cosPhi2);
			double dzp=-gc1.i*sin(gc1.zeta+t1*term1)+gc2.i*sin(gc2.zeta+t2*term2);
			*dp=2.0*dx*dxp+2.0*dy*dyp+2.0*dz*dzp;
		}
};

class ECubicFinder {
	public:
		static double collisionTime(CartCoords &c1,CartCoords &c2,EGCCoords &g1,EGCCoords &g2,double t1,double t2) {
			return 0.0;
		}
	private:
};

class ELinearFinder {
	public:
		static double collisionTime(CartCoords &c1,CartCoords &c2,EGCCoords &g1,EGCCoords &g2,double t1,double t2) {
			return 0.0;
		}
	private:
};

#endif

