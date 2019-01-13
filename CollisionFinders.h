// CollisionFinders.h
// This file contains various classes that can be used to find collisions
// between particles.

#include <stdio.h>
#include <math.h>
#include "GCPopulation.h"

#ifndef COLLISION_FINDERS
#define COLLISION_FINDERS

#ifdef STATS
class IterCounter {
	public:
		static void recordIters(int numIters,int outcome) {
			iterCnt[numIters][outcome]++;
		}
		static void recordCalls(int numCalls,int outcome) {
			callCnt[numCalls][outcome]++;
		}
		static void checkMade() {
			checkCnt++;
		}
		static void printNonZeroCounts(FILE *fout) {
			fprintf(fout,"total=%ld\n",checkCnt);
			fprintf(fout,"iters\n");
			for(int i=0; i<51; ++i) {
				long long int sum=sumRow(iterCnt[i]);
				if(sum!=0) fprintf(fout,"%d %lld\n",i,sum);
			}
			fprintf(fout,"calls\n");
			for(int i=0; i<1000; ++i) {
				long long int sum=sumRow(callCnt[i]);
				if(sum!=0) fprintf(fout,"%d %lld\n",i,sum);
			}
			fprintf(fout,"outcomes\n");
			for(int i=0; i<10; ++i) {
				long long int sum=sumColumn(iterCnt,i,50);
				if(sum!=0) fprintf(fout,"%d %lld\n",i,sum);
			}
		}
		static void printTables(FILE *fout) {
			for(int i=0; i<51; ++i) {
				if(sumRow(iterCnt[i])>0) {
					fprintf(fout,"-1 %d",i);
					for(int j=0; j<10; ++j) {
						fprintf(fout," %lld",iterCnt[i][j]);
					}
					fprintf(fout,"\n");
				}
			}
			for(int i=0; i<1000; ++i) {
				if(sumRow(callCnt[i])>0) {
					fprintf(fout,"-2 %d",i);
					for(int j=0; j<10; ++j) {
						fprintf(fout," %lld",callCnt[i][j]);
					}
					fprintf(fout,"\n");
				}
			}
		}
		static void clearCounts() {
			for(int i=0; i<51; ++i) {
				for(int j=0; j<10; ++j) {
					iterCnt[i][j]=0;
				}
			}
			for(int i=0; i<1000; ++i) {
				for(int j=0; j<10; ++j) {
					callCnt[i][j]=0;
				}
			}
			checkCnt=0;
		}

		static long long int sumRow(long long int row[]) {
			long long int sum=0;
			for(int i=0; i<10; ++i) {
				sum+=row[i];
			}
			return sum;
		}

		static long long int sumColumn(long long int table[][10],int col,int num) {
			long long int sum=0;
			for(int i=0; i<num; ++i) {
				sum+=table[i][col];
			}
			return sum;
		}

		static long long int iterCnt[][10];
		static long long int callCnt[][10];
		static long long int checkCnt;
};

long long int IterCounter::iterCnt[51][10];
long long int IterCounter::callCnt[1000][10];
long long int IterCounter::checkCnt;

#endif

class FullNewtonFinder {
	public:
		static double collisionTimeInitial(const CartCoords &cart1,const CartCoords &cart2,const GCCoords &gc1,const GCCoords &gc2,double radius1,double radius2,double dt) {
			return collisionTime(cart1,cart2,gc1,gc2,radius1,radius2,0.0,0.0,dt);
		}

		static double collisionTime(const CartCoords &cart1,const CartCoords &cart2,const GCCoords &gc1,const GCCoords &gc2,double radius1,double radius2,double time1,double time2,double dt) {
			double max_time=(time1>time2)?(time1):(time2);
			double dist_sqr,dist_sqr_prime;
			double rad_sqr=(radius1+radius2)*(radius1+radius2);
			GCDistancePrime(gc1,gc2,max_time-time1,max_time-time2,&dist_sqr,&dist_sqr_prime);
#ifdef STATS
			int calls=1;
			IterCounter::checkMade();
#endif
			if(dist_sqr<rad_sqr) {
				if(dist_sqr_prime<0.0) {
#ifdef STATS
					IterCounter::recordIters(0,4);
					IterCounter::recordCalls(calls,4);
#endif
					return(max_time+dt*0.0001);
				} else {
#ifdef STATS
					IterCounter::recordIters(0,2);
					IterCounter::recordCalls(calls,2);
#endif
					return(-2.0);
				}
			}
		
			double guess=max_time;
			if(dist_sqr_prime>0.0) {
#ifdef STATS
				IterCounter::recordIters(0,3);
				IterCounter::recordCalls(calls,3);
#endif
				return(-3.0);
			}
			for(int i=0; i<50; i++) {
				if(fabs(dist_sqr-rad_sqr)<1e-8*rad_sqr) {
					if((dist_sqr_prime>0.0) && (guess>0.0)) {
						while(dist_sqr_prime>0.0) {
							guess=0.5*(guess+max_time);
							GCDistancePrime(gc1,gc2,guess-time1,guess-time2,&dist_sqr,&dist_sqr_prime);
#ifdef STATS
							calls++;
#endif
						}
					} else {
		//				printf("Returning a guess of %e %e %e for %d %d\n",guess,dist_sqr,rad_sqr,p1,p2);
#ifdef STATS
						IterCounter::recordIters(i,0);
						IterCounter::recordCalls(calls,0);
#endif
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
#ifdef STATS
					calls++;
#endif
					cnt++;
				}
				if(cnt>=15) {
		//			printf("Cut in half too many times.\n");
#ifdef STATS
						IterCounter::recordIters(i,7);
						IterCounter::recordCalls(calls,7);
#endif
					return(-7.0);
				}
				guess=guess2;
				GCDistancePrime(gc1,gc2,guess-time1,guess-time2,&dist_sqr,&dist_sqr_prime);
#ifdef STATS
				calls++;
#endif
				if(guess<-5.0*dt) {
		//			printf("guess got too small %e\n",guess);
#ifdef STATS
						IterCounter::recordIters(i,5);
						IterCounter::recordCalls(calls,5);
#endif
					return(-5.0);
				}
				if(guess>5.0*dt) {
		//			printf("guess got too large %e\n",guess);
#ifdef STATS
						IterCounter::recordIters(i,6);
						IterCounter::recordCalls(calls,6);
#endif
					return(-6.0);
				}
			}
			printf("Didn't find a collision or jump out in 50 steps.  %e\n",guess);
			printf("Location 1 - %e %e %e %e %e %e\n",cart1.x,cart1.y,cart1.z,cart1.vx,cart1.vy,cart1.vz);
			printf("Location 2 - %e %e %e %e %e %e\n",cart2.x,cart2.y,cart2.z,cart2.vx,cart2.vy,cart2.vz);
#ifdef STATS
			IterCounter::recordIters(50,1);
			IterCounter::recordCalls(calls,1);
#endif
			return(-1.0);
		}
		static double GCDistance(const GCCoords &gc1,const GCCoords &gc2,double t1,double t2) {
			double dx=gc1.X-gc2.X-gc1.e*cos(gc1.phi+t1)+gc2.e*cos(gc2.phi+t2);
			double dy=(gc1.Y-1.5*gc1.X*t1)-(gc2.Y-1.5*gc2.X*t2)+2.0*gc1.e*sin(gc1.phi+t1)-2.0*gc2.e*sin(gc2.phi+t2);
			double dz=gc1.i*cos(gc1.zeta+t1)-gc2.i*cos(gc2.zeta+t2);
			return(dx*dx+dy*dy+dz*dz);
		}

		static void GCDistancePrime(const GCCoords &gc1,const GCCoords &gc2,double t1,double t2,double *d,double *dp) {
			double cosPhi1=cos(gc1.phi+t1);
			double cosPhi2=cos(gc2.phi+t2);
			double sinPhi1=sin(gc1.phi+t1);
			double sinPhi2=sin(gc2.phi+t2);
			double dx=gc1.X-gc2.X-gc1.e*cosPhi1+gc2.e*cosPhi2;
			double dy=(gc1.Y-1.5*gc1.X*t1)-(gc2.Y-1.5*gc2.X*t2)+2.0*gc1.e*sinPhi1-2.0*gc2.e*sinPhi2;
			double dz=gc1.i*cos(gc1.zeta+t1)-gc2.i*cos(gc2.zeta+t2);
			*d=dx*dx+dy*dy+dz*dz;
			double dxp=gc1.e*sinPhi1-gc2.e*sinPhi2;
			double dyp=-1.5*gc1.X+1.5*gc2.X+2.0*gc1.e*cosPhi1-2.0*gc2.e*cosPhi2;
			double dzp=-gc1.i*sin(gc1.zeta+t1)+gc2.i*sin(gc2.zeta+t2);
			*dp=2.0*dx*dxp+2.0*dy*dyp+2.0*dz*dzp;
		}
	private:
};

class FullNewtonFinder2 {
	public:
		static double collisionTime(const CartCoords &cart1,const CartCoords &cart2,const GCCoords &gc1,const GCCoords &gc2,double radius1,double radius2,double time1,double time2,double dt) {
			double max_time=(time1>time2)?(time1):(time2);
			double dist_sqr,dist_sqr_prime;
			double rad_sqr=(radius1+radius2)*(radius1+radius2);
			FullNewtonFinder::GCDistancePrime(gc1,gc2,max_time-time1,max_time-time2,&dist_sqr,&dist_sqr_prime);
#ifdef STATS
			int calls=1;
			IterCounter::checkMade();
#endif
			if(dist_sqr<rad_sqr) {
				if(dist_sqr_prime<0.0) {
#ifdef STATS
					IterCounter::recordIters(0,4);
					IterCounter::recordCalls(calls,4);
#endif
					return(max_time+dt*0.0001);
				} else {
#ifdef STATS
					IterCounter::recordIters(0,2);
					IterCounter::recordCalls(calls,2);
#endif
					return(-2.0);
				}
			}
		
			double guess=max_time;
			if(dist_sqr_prime>0.0) {
#ifdef STATS
				IterCounter::recordIters(0,3);
				IterCounter::recordCalls(calls,3);
#endif
				return(-3.0);
			}
			for(int i=0; i<50; i++) {
				if(fabs(dist_sqr-rad_sqr)<1e-8*rad_sqr) {
						//printf("Returning a guess of %e %e %e for %d %d\n",guess,dist_sqr,rad_sqr,p1,p2);
#ifdef STATS
					IterCounter::recordIters(i,0);
					IterCounter::recordCalls(calls,0);
#endif
					return(guess);
				}
				double step=-(dist_sqr-rad_sqr)/dist_sqr_prime;
				double guess2=guess+step;
				double dist_sqr2;
				FullNewtonFinder::GCDistancePrime(gc1,gc2,guess-time1,guess-time2,&dist_sqr2,&dist_sqr_prime);
#ifdef STATS
				calls++;
#endif
				int cnt=0;
				while((fabs(dist_sqr2-rad_sqr)>fabs(dist_sqr-rad_sqr) || dist_sqr_prime>0.0) && (cnt<15)) {
					//printf("%d %e %e %e %e\n",i,step,guess,dist_sqr2-rad_sqr,dist_sqr-rad_sqr);
					step*=0.5;
					guess2=guess+step;
					FullNewtonFinder::GCDistancePrime(gc1,gc2,guess2-time1,guess2-time2,&dist_sqr2,&dist_sqr_prime);
#ifdef STATS
					calls++;
#endif
					cnt++;
				}
				if(cnt>=15) {
					//printf("Cut in half too many times.\n");
#ifdef STATS
					IterCounter::recordIters(i,7);
					IterCounter::recordCalls(calls,7);
#endif
					return(-7.0);
				}
				guess=guess2;
				FullNewtonFinder::GCDistancePrime(gc1,gc2,guess-time1,guess-time2,&dist_sqr,&dist_sqr_prime);
#ifdef STATS
				calls++;
#endif
				if(guess<-5.0*dt) {
					//printf("guess got too small %e\n",guess);
#ifdef STATS
					IterCounter::recordIters(i,5);
					IterCounter::recordCalls(calls,5);
#endif
					return(-5.0);
				}
				if(guess>5.0*dt) {
					//printf("guess got too large %e\n",guess);
#ifdef STATS
					IterCounter::recordIters(i,6);
					IterCounter::recordCalls(calls,6);
#endif
					return(-6.0);
				}
			}
			printf("Didn't find a collision or jump out in 50 steps.  %e\n",guess);
			printf("Location 1 - %e %e %e %e %e %e\n",cart1.x,cart1.y,cart1.z,cart1.vx,cart1.vy,cart1.vz);
			printf("Location 2 - %e %e %e %e %e %e\n",cart2.x,cart2.y,cart2.z,cart2.vx,cart2.vy,cart2.vz);
#ifdef STATS
			IterCounter::recordIters(50,1);
			IterCounter::recordCalls(calls,1);
#endif
			return(-1.0);
		}
	private:
};

class LinearFinder {
	public:
		static double collisionTimeInitial(const CartCoords &cart1,const CartCoords &cart2,const GCCoords &gc1,const GCCoords &gc2,double radius1,double radius2,double dt) {
			double dx=cart1.x-cart2.x;
			double dy=cart1.y-cart2.y;
			double dz=cart1.z-cart2.z;
			double dvx=cart1.vx-cart2.vx;
			double dvy=cart1.vy-cart2.vy;
			double dvz=cart1.vz-cart2.vz;
			if(dx*dx+dy*dy+dz*dz<(radius1+radius2)*(radius1+radius2)*0.999 && 
				dx*dvx+dy*dvy+dz*dvz<0) return 1e-4*dt;
			double a=dvx*dvx+dvy*dvy+dvz*dvz;
			double b=2.0*(dx*dvx+dy*dvy+dz*dvz);
			double c=(dx*dx+dy*dy+dz*dz)-(radius1+radius2)*(radius1+radius2);
			double root=b*b-4.0*a*c;
			if(root<0) return -1;
			return (-b-sqrt(root))/(2.0*a);
		}

		static double collisionTime(const CartCoords &cc1,const CartCoords &cc2,const GCCoords &g1,const GCCoords &g2,double radius1,double radius2,double t1,double t2,double dt) {
			double max_time=(t1>t2)?(t1):(t2);
			GCCoords gc1(g1);
			GCCoords gc2(g2);
			gc1.advance(max_time-t1);
			gc2.advance(max_time-t2);
			CartCoords c1(gc1);
			CartCoords c2(gc2);
			double dx=c1.x-c2.x;
			double dy=c1.y-c2.y;
			double dz=c1.z-c2.z;
			double dvx=c1.vx-c2.vx;
			double dvy=c1.vy-c2.vy;
			double dvz=c1.vz-c2.vz;
			if(dx*dx+dy*dy+dz*dz<(radius1+radius2)*(radius1+radius2)*0.999 &&
				dx*dvx+dy*dvy+dz*dvz<0) return max_time+1e-4*dt;
			double a=dvx*dvx+dvy*dvy+dvz*dvz;
			double b=2.0*(dx*dvx+dy*dvy+dz*dvz);
			double c=(dx*dx+dy*dy+dz*dz)-(radius1+radius2)*(radius1+radius2);
			double root=b*b-4.0*a*c;
			if(root<0) return -1;
			double ret = (-b-sqrt(root))/(2.0*a);
			if(ret<max_time) return -1;
			return ret;
		}
		static double collisionOrNearestTime(CartCoords &c1,CartCoords &c2,double radius1,double radius2,int &moveMin) {
			double dx=c1.x-c2.x;
			double dy=c1.y-c2.y;
			double dz=c1.z-c2.z;
			double dvx=c1.vx-c2.vx;
			double dvy=c1.vy-c2.vy;
			double dvz=c1.vz-c2.vz;
			double a=dvx*dvx+dvy*dvy+dvz*dvz;
			double b=2.0*(dx*dvx+dy*dvy+dz*dvz);
			double c=(dx*dx+dy*dy+dz*dz)-(radius1+radius2)*(radius1+radius2);
			double root=b*b-4.0*a*c;
			if(root<0.0) {
				moveMin=1;
				return -b/(2.0*a);
			}
			moveMin=0;
			//printf("b=%e root=%e dqrt=%e a=%e c=%e\n",b,root,sqrt(root),a,c);
			return (-b-sqrt(root))/(2.0*a);
		}
	private:
};

class CubicFinder {
	public:
		static double collisionTime(CartCoords &cc1,CartCoords &cc2,GCCoords &g1,GCCoords &g2,double radius1,double radius2,double t1,double t2,double dt) {
#ifdef STATS
			IterCounter::checkMade();
#endif
			double max_time=(t1>t2)?(t1):(t2);
			GCCoords ig1(g1);
			GCCoords ig2(g2);
			ig1.advance(-t1);
			ig2.advance(-t2);
			CartCoords ic1(ig1);
			CartCoords ic2(ig2);
			double xi1[]={ic1.x,ic1.y,ic1.z};
			double xi2[]={ic2.x,ic2.y,ic2.z};
			double vi1[]={ic1.vx,ic1.vy,ic1.vz};
			double vi2[]={ic2.vx,ic2.vy,ic2.vz};
			GCCoords fg1(g1);
			GCCoords fg2(g2);
			fg1.advance(dt-t1);
			fg2.advance(dt-t2);
			CartCoords fc1(fg1);
			CartCoords fc2(fg2);
			double xf1[]={fc1.x,fc1.y,fc1.z};
			double xf2[]={fc2.x,fc2.y,fc2.z};
			double vf1[]={fc1.vx,fc1.vy,fc1.vz};
			double vf2[]={fc2.vx,fc2.vy,fc2.vz};
			double a1[3],b1[3],c1[3],d1[3];
			double a2[3],b2[3],c2[3],d2[3];
			cubicFit(xi1,vi1,xf1,vf1,a1,b1,c1,d1,dt);
			cubicFit(xi2,vi2,xf2,vf2,a2,b2,c2,d2,dt);
			double a[]={a1[0]-a2[0],a1[1]-a2[1],a1[2]-a2[2]};
			double b[]={b1[0]-b2[0],b1[1]-b2[1],b1[2]-b2[2]};
			double c[]={c1[0]-c2[0],c1[1]-c2[1],c1[2]-c2[2]};
			double d[]={d1[0]-d2[0],d1[1]-d2[1],d1[2]-d2[2]};
			double guess=max_time;
			double dist_sqr,dist_sqr_prime;
			double rad_sqr=(radius1+radius2)*(radius1+radius2);
			dist_sqr=cubicDistance(a,b,c,d,guess);
			dist_sqr_prime=cubicDistancePrime(a,b,c,d,guess);
			if(dist_sqr<rad_sqr) {
				if(dist_sqr_prime<0.0) {
#ifdef STATS
					IterCounter::recordIters(0,4);
#endif
					return(max_time+dt*0.0001);
				} else {
#ifdef STATS
					IterCounter::recordIters(0,2);
#endif
					return(-2.0);
				}
			}
		
			if(dist_sqr_prime>0.0) {
#ifdef STATS
				IterCounter::recordIters(0,3);
#endif
				return(-3.0);
			}
			for(int i=0; i<50; i++) {
				if(fabs(dist_sqr-rad_sqr)<1e-8*rad_sqr) {
						//printf("Returning a guess of %e %e %e for %d %d\n",guess,dist_sqr,rad_sqr,p1,p2);
#ifdef STATS
					IterCounter::recordIters(i,0);
#endif
					return(guess);
				}
				double step=-(dist_sqr-rad_sqr)/dist_sqr_prime;
				double guess2=guess+step;
				double dist_sqr2;
				dist_sqr2=cubicDistance(a,b,c,d,guess2);
				dist_sqr_prime=cubicDistancePrime(a,b,c,d,guess2);
				int cnt=0;
				while((fabs(dist_sqr2-rad_sqr)>fabs(dist_sqr-rad_sqr) || dist_sqr_prime>0.0) && (cnt<15)) {
					//printf("%d %e %e %e %e\n",i,step,guess,dist_sqr2-rad_sqr,dist_sqr-rad_sqr);
					step*=0.5;
					guess2=guess+step;
					dist_sqr2=cubicDistance(a,b,c,d,guess2);
					dist_sqr_prime=cubicDistancePrime(a,b,c,d,guess2);
					cnt++;
				}
				if(cnt>=15) {
					//printf("Cut in half too many times.\n");
#ifdef STATS
					IterCounter::recordIters(i,7);
#endif
					return(-7.0);
				}
				guess=guess2;
				dist_sqr=cubicDistance(a,b,c,d,guess);
				dist_sqr_prime=cubicDistancePrime(a,b,c,d,guess);
				if(guess<-5.0*dt) {
					//printf("guess got too small %e\n",guess);
#ifdef STATS
					IterCounter::recordIters(i,5);
#endif
					return(-5.0);
				}
				if(guess>5.0*dt) {
					//printf("guess got too large %e\n",guess);
#ifdef STATS
					IterCounter::recordIters(i,6);
#endif
					return(-6.0);
				}
			}
			printf("Didn't find a collision or jump out in 50 steps.  %e\n",guess);
			printf("Location 1 - %e %e %e %e %e %e\n",cc1.x,cc1.y,cc1.z,cc1.vx,cc1.vy,cc1.vz);
			printf("Location 2 - %e %e %e %e %e %e\n",cc2.x,cc2.y,cc2.z,cc2.vx,cc2.vy,cc2.vz);
#ifdef STATS
			IterCounter::recordIters(50,1);
#endif
			return -1.0;
		}

		static void cubicFit(double xi[],double vi[],double xf[],double vf[],double a[],double b[],double c[],double d[],double dt) {
			for(int i=0; i<3; ++i) {
				d[i]=xi[i];
				c[i]=vi[i];
				b[i]=(3*xf[i]-dt*vf[i]-2*c[i]*dt-3*d[i])/(dt*dt);
				a[i]=(xf[i]-b[i]*dt*dt-c[i]*dt-d[i])/(dt*dt*dt);
			}
		}

		static double cubicDistance(double a[],double b[],double c[],double d[],double t) {
			double dx=a[0]*t*t*t+b[0]*t*t+c[0]*t+d[0];
			double dy=a[1]*t*t*t+b[1]*t*t+c[1]*t+d[1];
			double dz=a[2]*t*t*t+b[2]*t*t+c[2]*t+d[2];
			return dx*dx+dy*dy+dz*dz;
		}

		static double cubicDistancePrime(double a[],double b[],double c[],double d[],double t) {
			double dx=2*(a[0]*t*t*t+b[0]*t*t+c[0]*t+d[0])*(3*a[0]*t*t+2*b[0]*t+c[0]);
			double dy=2*(a[1]*t*t*t+b[1]*t*t+c[1]*t+d[1])*(3*a[1]*t*t+2*b[1]*t+c[1]);
			double dz=2*(a[2]*t*t*t+b[2]*t*t+c[2]*t+d[2])*(3*a[2]*t*t+2*b[2]*t+c[2]);
			return dx+dy+dz;
		}
	private:
};

template<class GCType>
class FullLinearFinder {
	public:
		static double collisionTime(CartCoords &c1,CartCoords &c2,GCType &g1,GCType &g2,double radius1,double radius2,double t1,double t2,double dt) {
#ifdef STATS
			int calls=1;
			IterCounter::checkMade();
#endif
			double max_time=(t1>t2)?(t1):(t2);
			double rad_sqr=(radius1+radius2)*(radius1+radius2);
			GCType gc1(g1);
			GCType gc2(g2);
			gc1.advance(max_time-t1);
			gc2.advance(max_time-t2);
			CartCoords cc1(gc1);
			CartCoords cc2(gc2);
			double dist_sqr,dot;
			dist_sqr=cartDistance(cc1,cc2,dot);
			if(dist_sqr<=rad_sqr*1.00000001) {
				if(dot<0.0) {
#ifdef STATS
					IterCounter::recordIters(0,4);
					IterCounter::recordCalls(calls,4);
#endif
					if(max_time==0.0) {
						return(max_time+dt*0.0001);
					} else {
						double dvx=cc1.vx-cc2.vx;
						double dvy=cc1.vy-cc2.vy;
						double dvz=cc1.vz-cc2.vz;
						double vel=sqrt(dvx*dvx+dvy*dvy+dvz*dvz);
						double minRad=(radius1<radius2)?radius1:radius2;
						//printf("%e %e\n",dt*0.0001,0.001*minRad/vel);
						return(max_time+0.001*minRad/vel);
					}
				} else {
#ifdef STATS
					IterCounter::recordIters(0,2);
					IterCounter::recordCalls(calls,2);
#endif
					return -2;
				}
			}
			if(dot>0.0) {
#ifdef STATS
				IterCounter::recordIters(0,3);
				IterCounter::recordCalls(calls,3);
#endif
				return(-3.0);
			}
			return collisionTimeLoop(cc1,cc2,gc1,gc2,radius1,radius2,dt,rad_sqr,max_time);
		}

		static double collisionTimeInitial(CartCoords &c1,CartCoords &c2,GCType &g1,GCType &g2,double radius1,double radius2,double dt) {
#ifdef STATS
			int calls=1;
			IterCounter::checkMade();
#endif
			double rad_sqr=(radius1+radius2)*(radius1+radius2);
			GCType gc1(g1);
			GCType gc2(g2);
			CartCoords cc1(c1);
			CartCoords cc2(c2);
			double dist_sqr,dot;
			dist_sqr=cartDistance(cc1,cc2,dot);
			if(dist_sqr<=rad_sqr*1.00000001) {
				if(dot<0.0) {
#ifdef STATS
					IterCounter::recordIters(0,4);
					IterCounter::recordCalls(calls,4);
#endif
					return(dt*0.0001);
				} else {
#ifdef STATS
					IterCounter::recordIters(0,2);
					IterCounter::recordCalls(calls,2);
#endif
					return -2;
				}
			}
			if(dot>0.0) {
#ifdef STATS
				IterCounter::recordIters(0,3);
				IterCounter::recordCalls(calls,3);
#endif
				return(-3.0);
			}
			return collisionTimeLoop(cc1,cc2,gc1,gc2,radius1,radius2,dt,rad_sqr,0);
		}

		static double collisionTimeLoop(CartCoords &cc1,CartCoords &cc2,GCType &gc1,GCType &gc2,double radius1,double radius2,double dt,double rad_sqr,double max_time) {
			double guess=max_time;
			int minCount=0;
			double guess2=dt;
			int moveMin=0;
			double dot=0;
#ifdef STATS
			calls++;
#endif
			for(int i=0; i<50; ++i) {
				double dist_sqr=cartDistance(cc1,cc2,dot);
				if(fabs(dist_sqr-rad_sqr)<1e-8*rad_sqr && dot<0.0) {
#ifdef STATS
					IterCounter::recordIters(i,0);
					IterCounter::recordCalls(calls,0);
#endif
					return(guess);
				}
				if(moveMin && dist_sqr>2*rad_sqr) {
#ifdef STATS
					IterCounter::recordIters(i,8);
					IterCounter::recordCalls(calls,8);
#endif
					return -8.0;
				}
				guess2=LinearFinder::collisionOrNearestTime(cc1,cc2,radius1,radius2,moveMin);
				minCount+=moveMin;
				if(minCount>25) {
#ifdef STATS
					IterCounter::recordIters(i,9);
					IterCounter::recordCalls(calls,9);
#endif
					return -9.0;
				}
				if(fabs(guess2)<1e-8*dt && moveMin) {
					//printf("Not moving\n");
#ifdef STATS
					IterCounter::recordIters(i,8);
					IterCounter::recordCalls(calls,8);
#endif
					return -8.0;
				}
				gc1.advance(guess2);
				gc2.advance(guess2);
				guess+=guess2;
				cc1.set(gc1);
				cc2.set(gc2);
#ifdef STATS
				calls++;
#endif
				double dist_sqr2=cartDistance(cc1,cc2,dot);
				int cnt=0;
				while((fabs(dist_sqr2-rad_sqr)>fabs(dist_sqr-rad_sqr) || dot>0.0) && cnt<15) {
					//printf("check %e %e %e %d %e %e\n",dist_sqr2,dist_sqr,dot,cnt,guess,guess2);
					guess2*=0.5;
					guess-=guess2;
					gc1.advance(-guess2);
					gc2.advance(-guess2);
					cc1.set(gc1);
					cc2.set(gc2);
#ifdef STATS
					calls++;
#endif
					dist_sqr2=cartDistance(cc1,cc2,dot);
					cnt++;
				}
				if(cnt>=15) {
					//printf("Not moving\n");
#ifdef STATS
					IterCounter::recordIters(i,7);
					IterCounter::recordCalls(calls,7);
#endif
					return -7.0;
				}
				if(guess<-5.0*dt) {
					//printf("guess got too small %e\n",guess);
#ifdef STATS
					IterCounter::recordIters(i,5);
					IterCounter::recordCalls(calls,5);
#endif
					return -5.0;
				}
				if(guess>5.0*dt) {
					//printf("guess got too large %e\n",guess);
#ifdef STATS
					IterCounter::recordIters(i,6);
					IterCounter::recordCalls(calls,6);
#endif
					return -6.0;
				}
			}
			printf("Didn't find a collision or jump out in 50 steps.  %e\n",guess);
			printf("Location 1 - %e %e %e %e %e %e\n",cc1.x,cc1.y,cc1.z,cc1.vx,cc1.vy,cc1.vz);
			printf("Location 2 - %e %e %e %e %e %e\n",cc2.x,cc2.y,cc2.z,cc2.vx,cc2.vy,cc2.vz);
#ifdef STATS
			IterCounter::recordIters(50,1);
			IterCounter::recordCalls(calls,1);
#endif
			return -1.0;
		}

		static double cartDistance(CartCoords &cc1,CartCoords &cc2,double &dot) {
			double dx=cc1.x-cc2.x;
			double dy=cc1.y-cc2.y;
			double dz=cc1.z-cc2.z;
			double dvx=cc1.vx-cc2.vx;
			double dvy=cc1.vy-cc2.vy;
			double dvz=cc1.vz-cc2.vz;
			dot=dx*dvx+dy*dvy+dz*dvz;
			return dx*dx+dy*dy+dz*dz;
		}

	private:
};

class CubicLinearFinder {
	public:
		static double collisionTime(CartCoords &cart1,CartCoords &cart2,GCCoords &g1,GCCoords &g2,double radius1,double radius2,double t1,double t2,double dt) {
#ifdef STATS
			IterCounter::checkMade();
#endif
			double max_time=(t1>t2)?(t1):(t2);
			double rad_sqr=(radius1+radius2)*(radius1+radius2);
			GCCoords ig1(g1);
			GCCoords ig2(g2);
			ig1.advance(max_time-t1);
			ig2.advance(max_time-t2);
			CartCoords ic1(ig1);
			CartCoords ic2(ig2);
			double dx=ic1.x-ic2.x;
			double dy=ic1.y-ic2.y;
			double dz=ic1.z-ic2.z;
			double dvx=ic1.vx-ic2.vx;
			double dvy=ic1.vy-ic2.vy;
			double dvz=ic1.vz-ic2.vz;
			if(dx*dx+dy*dy+dz*dz<=rad_sqr*1.00000001) {
				if(dx*dvx+dy*dvy+dz*dvz<0.0) {
#ifdef STATS
					IterCounter::recordIters(0,4);
#endif
					return(max_time+dt*0.0001);
				} else {
#ifdef STATS
					IterCounter::recordIters(0,2);
#endif
					return -2;
				}
			}
			if(dx*dvx+dy*dvy+dz*dvz>0.0) {
#ifdef STATS
				IterCounter::recordIters(0,3);
#endif
				return(-3.0);
			}
			ig1.advance(-max_time);
			ig2.advance(-max_time);
			ic1.set(ig1);
			ic2.set(ig2);
			double xi1[]={ic1.x,ic1.y,ic1.z};
			double xi2[]={ic2.x,ic2.y,ic2.z};
			double vi1[]={ic1.vx,ic1.vy,ic1.vz};
			double vi2[]={ic2.vx,ic2.vy,ic2.vz};
			GCCoords fg1(g1);
			GCCoords fg2(g2);
			fg1.advance(dt-t1);
			fg2.advance(dt-t2);
			CartCoords fc1(fg1);
			CartCoords fc2(fg2);
			double xf1[]={fc1.x,fc1.y,fc1.z};
			double xf2[]={fc2.x,fc2.y,fc2.z};
			double vf1[]={fc1.vx,fc1.vy,fc1.vz};
			double vf2[]={fc2.vx,fc2.vy,fc2.vz};
			double a1[3],b1[3],c1[3],d1[3];
			double a2[3],b2[3],c2[3],d2[3];
			CubicFinder::cubicFit(xi1,vi1,xf1,vf1,a1,b1,c1,d1,dt);
			CubicFinder::cubicFit(xi2,vi2,xf2,vf2,a2,b2,c2,d2,dt);
			double a[]={a1[0]-a2[0],a1[1]-a2[1],a1[2]-a2[2]};
			double b[]={b1[0]-b2[0],b1[1]-b2[1],b1[2]-b2[2]};
			double c[]={c1[0]-c2[0],c1[1]-c2[1],c1[2]-c2[2]};
			double d[]={d1[0]-d2[0],d1[1]-d2[1],d1[2]-d2[2]};
			double guess=max_time;
			CartCoords cc1,cc2;
			int minCount=0;
			for(int i=0; i<50; ++i) {
				double dist_sqr=CubicFinder::cubicDistance(a,b,c,d,guess);
				double dist_sqr_prime=CubicFinder::cubicDistancePrime(a,b,c,d,guess);
				if(fabs(dist_sqr-rad_sqr)<1e-8*rad_sqr && dist_sqr_prime<0.0) {
#ifdef STATS
					IterCounter::recordIters(i,0);
#endif
					return(guess);
				}
				setCartFromCubic(a1,b1,c1,d1,guess,cc1);
				setCartFromCubic(a2,b2,c2,d2,guess,cc2);
				int moveMin;
				double guess2=LinearFinder::collisionOrNearestTime(cc1,cc2,radius1,radius2,moveMin);
				minCount+=moveMin;
				if(minCount>25) {
#ifdef STATS
					IterCounter::recordIters(i,9);
#endif
					return -9.0;
				}
				if(fabs(guess2)<1e-8*dt && moveMin) {
					//printf("Not moving\n");
#ifdef STATS
					IterCounter::recordIters(i,8);
#endif
					return -8.0;
				}
				double dist_sqr2=CubicFinder::cubicDistance(a,b,c,d,guess+guess2);
				dist_sqr_prime=CubicFinder::cubicDistancePrime(a,b,c,d,guess+guess2);
				int cnt=0;
				while((fabs(dist_sqr2-rad_sqr)>fabs(dist_sqr-rad_sqr) || dist_sqr_prime>0.0) && cnt<15) {
					//printf("check %e %e %e %d %e %e\n",dist_sqr2,dist_sqr,dot,cnt,guess,guess2);
					guess2*=0.5;
					dist_sqr2=CubicFinder::cubicDistance(a,b,c,d,guess+guess2);
					dist_sqr_prime=CubicFinder::cubicDistancePrime(a,b,c,d,guess+guess2);
					cnt++;
				}
				if(cnt>=15) {
					//printf("Not moving\n");
#ifdef STATS
					IterCounter::recordIters(i,7);
#endif
					return -7.0;
				}
				guess+=guess2;
				if(guess<-5.0*dt) {
					//printf("guess got too small %e\n",guess);
#ifdef STATS
					IterCounter::recordIters(i,5);
#endif
					return(-5.0);
				}
				if(guess>5.0*dt) {
					//printf("guess got too large %e\n",guess);
#ifdef STATS
					IterCounter::recordIters(i,6);
#endif
					return(-6.0);
				}
			}
			printf("Didn't find a collision or jump out in 50 steps.  %e\n",guess);
			printf("Location 1 - %e %e %e %e %e %e\n",cart1.x,cart1.y,cart1.z,cart1.vx,cart1.vy,cart1.vz);
			printf("Location 2 - %e %e %e %e %e %e\n",cart2.x,cart2.y,cart2.z,cart2.vx,cart2.vy,cart2.vz);
#ifdef STATS
			IterCounter::recordIters(50,1);
#endif
			return -1.0;
		}

		static void setCartFromCubic(double a[],double b[],double c[],double d[],double t,CartCoords &cc) {
			cc.x=a[0]*t*t*t+b[0]*t*t+c[0]*t+d[0];
			cc.y=a[1]*t*t*t+b[1]*t*t+c[1]*t+d[1];
			cc.z=a[2]*t*t*t+b[2]*t*t+c[2]*t+d[2];
			cc.vx=3*a[0]*t*t+2*b[0]*t+c[0];
			cc.vy=3*a[1]*t*t+2*b[1]*t+c[1];
			cc.vz=3*a[2]*t*t+2*b[2]*t+c[2];
		}

	private:
};

#endif

