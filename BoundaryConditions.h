// BoundaryConditions.h
// This file contains some classes for different boundary conditions.  These
// classes are one of the main places where code for multiple processors shows
// up, because without gravity using boundaing layers the collisions run like
// they are on a single machine for each timestep.  At the end of the timestep
// the talk to one another about which particles have moved across the lines.
// A significant question here is whether I should make separate versions for
// the sinlge processor BCs or if that should be special code in a more
// complete version.

#ifndef BOUNDARY_CONDITIONS
#define BOUNDARY_CONDITIONS

#include <vector>
#include <tuple>
//#include "GCPopulation.h"
#include "ProcessorCommunication.h"
#include "Coordinates.h"

using std::tuple;

class OpenBounds {
	public:
		template<class Population>
		void apply(Population &pop) {}

		bool checkIfIn(double x,double y) {
			return true;
		}

		double getArea() {
			return 1e100;
		}

		double getMinX() { return -1e100; }
		double getMaxX() { return 1e100; }
		double getMinY() { return -1e100; }
		double getMaxY() { return 1e100; }
		double getMinZ() { return -1e100; }
		double getMaxZ() { return 1e100; }
};

class ShearedPeriodic {
	public:
		ShearedPeriodic(double sx,double sy,double sz,double sh):sizex(sx),sizey(sy),sizez(sz),shear(sh) {}

		template<class Population>
		void apply(Population &pop) {
			int nb = pop.getNumBodies();
			#pragma omp parallel for schedule(static)
			for(int i=0; i<nb; ++i) {
				while(pop.getx(i)<-sizex) {
					pop.setx(i,pop.getx(i)+(2.0*sizex));
					pop.setvy(i,pop.getvx(i)+shear);
				} 
				while(pop.getx(i)>sizex) {
					pop.setx(i,pop.getx(i)-(2.0*sizex));
					pop.setvy(i,pop.getvx(i)-shear);
				} 
				while(pop.gety(i)<-sizey) {
					pop.sety(i,pop.gety(i)+(2.0*sizey));
				} 
				while(pop.gety(i)>sizey) {
					pop.sety(i,pop.gety(i)-(2.0*sizey));
				}
				while(pop.getz(i)<-sizez) {
					pop.setz(i,pop.getz(i)+(2.0*sizez));
				} 
				while(pop.getz(i)>sizez) {
					pop.setz(i,pop.getz(i)-(2.0*sizez));
				} 
			}
		}

		bool checkIfIn(double x,double y) {
			return (x>=-sizex && x<=sizex && y>=-sizey && y<=sizey);
		}

		double getArea() {
			return (2.0*sizex)*(2.0*sizey);
		}

		double getMinX() { return -sizex; }
		double getMaxX() { return sizex; }
		double getMinY() { return -sizey; }
		double getMaxY() { return sizey; }
		double getMinZ() { return -sizez; }
		double getMaxZ() { return sizez; }
	private:
		double sizex,sizey,sizez;
		double shear;
};

class SlidingBrick {
	public:
		SlidingBrick(double minX,double maxX,double minY,double maxY):minx(minX),maxx(maxX),miny(minY),maxy(maxY) {}

		template<class Population>
		void apply(Population &pop) {
			int nb = pop.getNumBodies();
			#pragma omp parallel for schedule(static)
			for(int ii=0; ii<nb; ++ii) {
				ParticleIndex pi = {ii};
				while(pop.getx(pi)<minx) {
					pop.setx(pi,pop.getx(pi)+(maxx-minx));
					pop.setvy(pi,pop.getvy(pi)-2.0*GCCoords::A0*(maxx-minx));
					pop.adjustAfterForce(pi);
				} 
				while(pop.getx(pi)>maxx) {
					pop.setx(pi,pop.getx(pi)-(maxx-minx));
					pop.setvy(pi,pop.getvy(pi)+2.0*GCCoords::A0*(maxx-minx));
					pop.adjustAfterForce(pi);
				}
				while(pop.gety(pi)<miny) {
					pop.sety(pi,pop.gety(pi)+(maxy-miny));
					pop.adjustAfterForce(pi);
				} 
				while(pop.gety(pi)>maxy) {
					pop.sety(pi,pop.gety(pi)-(maxy-miny));
					pop.adjustAfterForce(pi);
				}
			}
		}

		bool checkIfIn(double x,double y) {
			return (x>=minx && x<=maxx && y>=miny && y<=maxy);
		}

		double getArea() {
			return (maxx-minx)*(maxy-miny);
		}

		template<class Population>
		void offset(Population &pop, ParticleIndex pi, double offsetX, double offsetY, CartCoords &cc) {
			cc = pop.getCart()[pi.i];
			cc.y += offsetY;
			if (offsetX != 0.0) {
				cc.vy -= 2.0*GCCoords::A0*offsetX;
				cc.x += offsetX;
			}
		}

		double getMinX() { return minx; }
		double getMaxX() { return maxx; }
		double getMinY() { return miny; }
		double getMaxY() { return maxy; }
	private:
		double minx,maxx,miny,maxy;
};

// This is my normal boundary for a perturbed system.  It is quasi-periodic
// in azimuth with a systematic adjustment to the phase angle applied.
template<class GCType>
class PeriodicWithPhiShift {
	public:
		PeriodicWithPhiShift(double minX,double maxX,double minY,double maxY):minx(minX),maxx(maxX),miny(minY),maxy(maxY) {
			if(minx > 0.0) edgex = minx;
			else edgex = maxx;
		}

		template<class Population>
		void apply(Population &pop) {
			miny+=GCType::Ydot(edgex)*pop.getTimeStep();
			maxy+=GCType::Ydot(edgex)*pop.getTimeStep();
			int nb = pop.getNumBodies();
			#pragma omp parallel for schedule(static)
			for(int ii=0; ii<nb; ++ii) {
				ParticleIndex pi = {ii};
				while(pop.gety(pi)<miny) {
					pop.setY(pi,pop.getY(pi)+(maxy-miny));
					pop.setPhi(pi,pop.getPhi(pi)+GCType::phidot(pop.getX(pi))*(maxy-miny)/GCType::Ydot(pop.getX(pi)));
					pop.setCartAfterForce(pi);
				}
				while(pop.gety(pi)>maxy) {
					pop.setY(pi,pop.getY(pi)-(maxy-miny));
					pop.setPhi(pi,pop.getPhi(pi)-GCType::phidot(pop.getX(pi))*(maxy-miny)/GCType::Ydot(pop.getX(pi)));
					pop.setCartAfterForce(pi);
				}
			}
		}

		bool checkIfIn(double x,double y) {
			return (x>=minx && x<=maxx && y>=miny && y<=maxy);
		}

		double getArea() {
			return (maxx-minx)*(maxy-miny);
		}

		template<class Population>
		void offset(Population &pop, ParticleIndex pi, double offsetX, double offsetY, CartCoords &cc) {
			GCType gc = pop.getGC()[pi.i];
			gc.Y += offsetY;
			gc.phi += GCType::phidot(pop.getX(pi))*offsetY/GCType::Ydot(pop.getX(pi));
			cc.set(gc);
		}

		double getMinX() { return minx; }
		double getMaxX() { return maxx; }
		double getMinY() { return miny; }
		double getMaxY() { return maxy; }
	private:
		double minx,maxx,miny,maxy;
		double edgex;
};

// This BC keeps an azimuthal extent that is a single orbital period around.
// That means that is isn't square as the length gets smaller as one moves
// closer to the perturber.  This should probably also do some strange things
// to handle the case of when the cell is moving past the perturber.  The
// populations at the two edges won't match well during that time so I will
// want the boundaries "off" at that time.
class SingleOrbitAzimuthal {
	public:
		SingleOrbitAzimuthal(double minX,double maxX,double minY):minx(minX),maxx(maxX),miny(minY) {}

		template<class Population>
		void apply(Population &pop) {
			miny-=2.0*GCCoords::A0*minx*pop.getTimeStep();
			double deltay=4.0*3.141592654*GCCoords::A0*maxx;
			double tmpy=miny;
			while(tmpy<-3.141592654) tmpy+=6.283185307;
			if(tmpy<-2.0*deltay || tmpy>deltay) {
				int nb = pop.getNumBodies();
				#pragma omp parallel for schedule(static)
				for(int ii=0; ii<nb; ++ii) {
					ParticleIndex pi = {ii};
					while(pop.gety(pi)<miny) {
						pop.setY(pi,pop.getY(pi)+4.0*3.141592654*GCCoords::A0*pop.getX(pi));
						pop.setCartAfterForce(pi);
					}
					while(pop.gety(pi)>miny+4.0*3.141592654*GCCoords::A0*pop.getx(pi)) {
						pop.setY(pi,pop.getY(pi)-4.0*3.141592654*GCCoords::A0*pop.getX(pi));
						pop.setCartAfterForce(pi);
					}
				}
			}
		}

		bool checkIfIn(double x,double y) {
			return (x>=minx && x<=maxx && y>=miny && y<=miny+4.0*3.141592654*GCCoords::A0*x);
		}

		double getArea() {
			return (getMaxX()-getMinX())*(getMaxY()-getMinY());
//			return (getMaxX()-getMinX())*(getMaxY()-getMinY())-0.5*(getMaxX()-getMinX())*(4.0*3.141592654*GCCoords::A0*(getMaxX()-getMinX()));
		}

		double getMinX() { return minx; }
		double getMaxX() { return maxx; }
		double getMinY() { return miny; }
		double getMaxY() { return miny+4.0*3.141592654*GCCoords::A0*maxx; }
	private:
		double minx,maxx;
		// I don't store maxy because maxy(x)=miny+4*PI*A0*x
		double miny;
};

// This is my normal boundary for a perturbed system.  It is quasi-periodic
// in azimuth with a systematic adjustment to the phase angle applied.
class FixedPeriodic {
	public:
		FixedPeriodic(double minX,double maxX,double minY,double maxY,bool zero=false)
			:minx(minX),maxx(maxX),miny(minY),maxy(maxY),zeroOnWrap(zero) {}

		template<class Population>
		void apply(Population &pop) {
			int nb = pop.getNumBodies();
			#pragma omp parallel for schedule(static)
			for(int ii=0; ii<nb; ii++) {
				ParticleIndex pi = {ii};
				while(pop.gety(pi)<miny) {
					pop.setY(pi,pop.getY(pi)+(maxy-miny));
					if(zeroOnWrap) pop.sete(pi,0);
					pop.setCartAfterForce(pi);
				}
				while(pop.gety(pi)>maxy) {
					pop.setY(pi,pop.getY(pi)-(maxy-miny));
					pop.setCartAfterForce(pi);
				}
			}
		}

		bool checkIfIn(double x,double y) {
			return (x>=minx && x<=maxx && y>=miny && y<=maxy);
		}

		double getArea() {
			return (maxx-minx)*(maxy-miny);
		}

		double getMinX() { return minx; }
		double getMaxX() { return maxx; }
		double getMinY() { return miny; }
		double getMaxY() { return maxy; }
	private:
		double minx,maxx,miny,maxy;
		bool zeroOnWrap;
};

#ifdef PARALLEL

// Same as earlier, but now written to handle running in parallel.
template<class SpacingMeasure>
class ParallelPeriodicWithPhiShift {
	public:
		ParallelPeriodicWithPhiShift(double minX,double maxX,double minY,double maxY,ProcessorCommunication &procComm,SpacingMeasure &sm):minxTot(minX),maxxTot(maxX),minyTot(minY),maxyTot(maxY),pc(procComm),spaceMeasure(sm) {
			minx=minxTot;
			maxx=maxxTot;
			double size=(maxyTot-minyTot)/pc.getNumProcesses();
			miny=minyTot+pc.getProcessNum()*size;
			maxy=minyTot+(pc.getProcessNum()+1)*size;
		}

		template<class Population>
		void apply(Population &pop) {
			double bufferSize=3.0*spaceMeasure.getMinSpacing();
			miny-=2.0*GCCoords::A0*minx*pop.getTimeStep();
			maxy-=2.0*GCCoords::A0*minx*pop.getTimeStep();
			minyTot-=2.0*GCCoords::A0*minx*pop.getTimeStep();
			maxyTot-=2.0*GCCoords::A0*minx*pop.getTimeStep();
			int numMoved[2]={0,0};
			int numBuffer[2]={0,0};
			commBuffer[0].resize(2);
			commBuffer[1].resize(2);
			removed.resize(0);
			sizeBuffer.resize(2);

			double deltay=4.0*3.141592654*GCCoords::A0*maxx;
			double tmpy=miny;
			while(tmpy<-3.141592654) tmpy+=6.283185307;
			bool wrapAround=(tmpy<-2.0*deltay || tmpy>deltay);

//			bool wrapAround=true;

			printf("%d - Applying BC %e %e %e %d\n",pc.getProcessNum(),miny,maxy,bufferSize,pop.getNumReal());

			// Only go up to numReal.  Others are thrown away.
			for(int i=0; i<pop.getNumReal(); ++i) {
				// Check downstream.
//				if(i%10000==0) printf("i[%d]=%e\n",i,pop.gety(i));
				if(pop.gety(i)<miny) {
					if(pc.getProcessNum()==0 && wrapAround) {
						pop.setY(i,pop.getY(i)+(maxyTot-minyTot));
						pop.setPhi(i,pop.getPhi(i)-(maxyTot-minyTot)/(2.0*GCCoords::A0*pop.getX(i)));
						pop.setCartAfterForce(i);
					}
					if(wrapAround || pc.getProcessNum()>0) {
						commBuffer[0].resize(commBuffer[0].size()+7);
						if(numBuffer[0]>0) {
							for(int j=0; j<7; ++j) {
								commBuffer[0][2+7*numMoved[0]+7*numBuffer[0]+j]=
								commBuffer[0][2+7*numMoved[0]+j];
							}
						}
						commBuffer[0][2+7*numMoved[0]]=pop.getx(i);
						commBuffer[0][3+7*numMoved[0]]=pop.gety(i);
						commBuffer[0][4+7*numMoved[0]]=pop.getz(i);
						commBuffer[0][5+7*numMoved[0]]=pop.getvx(i);
						commBuffer[0][6+7*numMoved[0]]=pop.getvy(i);
						commBuffer[0][7+7*numMoved[0]]=pop.getvz(i);
						commBuffer[0][8+7*numMoved[0]]=pop.getRadius(i);
						numMoved[0]++;
						removed.push_back(i);
					}
				} else if(pop.gety(i)<miny+bufferSize && pc.getProcessNum()>0) {
					commBuffer[0].resize(commBuffer[0].size()+7);
					commBuffer[0][2+7*numMoved[0]+7*numBuffer[0]]=pop.getx(i);
					commBuffer[0][3+7*numMoved[0]+7*numBuffer[0]]=pop.gety(i);
					commBuffer[0][4+7*numMoved[0]+7*numBuffer[0]]=pop.getz(i);
					commBuffer[0][5+7*numMoved[0]+7*numBuffer[0]]=pop.getvx(i);
					commBuffer[0][6+7*numMoved[0]+7*numBuffer[0]]=pop.getvy(i);
					commBuffer[0][7+7*numMoved[0]+7*numBuffer[0]]=pop.getvz(i);
					commBuffer[0][8+7*numMoved[0]+7*numBuffer[0]]=pop.getRadius(i);
					numBuffer[0]++;
				} else if(pop.gety(i)>maxy) {
				// Check upstream.
					if(pc.getProcessNum()==pc.getNumProcesses()-1 && wrapAround) {
						pop.setY(i,pop.getY(i)-(maxyTot-minyTot));
						pop.setPhi(i,pop.getPhi(i)+(maxyTot-minyTot)/(2.0*GCCoords::A0*pop.getX(i)));
						pop.setCartAfterForce(i);
					}
					if(wrapAround || pc.getProcessNum()<pc.getNumProcesses()-1) {
						commBuffer[1].resize(commBuffer[1].size()+7);
						if(numBuffer[1]>0) {
							for(int j=0; j<7; ++j) {
								commBuffer[1][2+7*numMoved[1]+7*numBuffer[1]+j]=
								commBuffer[1][2+7*numMoved[1]+j];
							}
						}
						commBuffer[1][2+7*numMoved[1]]=pop.getx(i);
						commBuffer[1][3+7*numMoved[1]]=pop.gety(i);
						commBuffer[1][4+7*numMoved[1]]=pop.getz(i);
						commBuffer[1][5+7*numMoved[1]]=pop.getvx(i);
						commBuffer[1][6+7*numMoved[1]]=pop.getvy(i);
						commBuffer[1][7+7*numMoved[1]]=pop.getvz(i);
						commBuffer[1][8+7*numMoved[1]]=pop.getRadius(i);
						numMoved[1]++;
						removed.push_back(i);
					}
				} else if(pop.gety(i)>maxy-bufferSize && pc.getProcessNum()<pc.getNumProcesses()-1) {
					commBuffer[1].resize(commBuffer[1].size()+7);
					commBuffer[1][2+7*numMoved[1]+7*numBuffer[1]]=pop.getx(i);
					commBuffer[1][3+7*numMoved[1]+7*numBuffer[1]]=pop.gety(i);
					commBuffer[1][4+7*numMoved[1]+7*numBuffer[1]]=pop.getz(i);
					commBuffer[1][5+7*numMoved[1]+7*numBuffer[1]]=pop.getvx(i);
					commBuffer[1][6+7*numMoved[1]+7*numBuffer[1]]=pop.getvy(i);
					commBuffer[1][7+7*numMoved[1]+7*numBuffer[1]]=pop.getvz(i);
					commBuffer[1][8+7*numMoved[1]+7*numBuffer[1]]=pop.getRadius(i);
					numBuffer[1]++;
				}

				// Check inward if neighbors that direction.
				// !!! Not written yet.

				// Check outward if neighbors that direction.
				// !!! Not written yet.
			}
			// Send data to neighbors.
			printf("%d - Sending in BC %d %d %d %d\n",pc.getProcessNum(),numMoved[0],numBuffer[0],numMoved[1],numBuffer[1]);
//			printf("Check 0 %d - %d\n",commBuffer[0].size(),2+7*(numMoved[0]+numBuffer[0]));
//			printf("Check 1 %d - %d\n",commBuffer[1].size(),2+7*(numMoved[1]+numBuffer[1]));
//			if(pc.getProcessNum()==1) {
//				printf("Sending 0 - ");
//				for(int i=0; i<commBuffer[0].size(); i++) printf("%e ",commBuffer[0][i]);
//				printf("\n");
//				printf("Sending 1 - ");
//				for(int i=0; i<2+7*(numMoved[1]+numBuffer[1]); i++) printf("%e ",commBuffer[1][i]);
//				printf("\n");
//			}
			commBuffer[0][0]=numMoved[0];
			commBuffer[0][1]=numBuffer[0];
			commBuffer[1][0]=numMoved[1];
			commBuffer[1][1]=numBuffer[1];

			// Read data from neighbors.
//			printf("%d - Reading in BC\n",pc.getProcessNum());
			int numRead[2];
			int error[2];

			// Send right, read left.
			sizeBuffer[0]=numMoved[0];
			sizeBuffer[1]=numBuffer[0];
			pc.sendToNeighbor(3,sizeBuffer);
			pc.sendToNeighbor(3,commBuffer[0]);
			error[0]=pc.readFromNeighbor(5,sizeBuffer,numRead[0]);
			unsigned int size=2+7*(int)(sizeBuffer[0]+sizeBuffer[1]);
			if(size>commBuffer[0].size())
				commBuffer[0].resize(size);
			error[0]=pc.readFromNeighbor(5,commBuffer[0],numRead[0]);

			// Send left, read right.
			sizeBuffer[0]=numMoved[1];
			sizeBuffer[1]=numBuffer[1];
			pc.sendToNeighbor(5,sizeBuffer);
			pc.sendToNeighbor(5,commBuffer[1]);
			error[1]=pc.readFromNeighbor(3,sizeBuffer,numRead[1]);
			size=2+7*(int)(sizeBuffer[0]+sizeBuffer[1]);
			if(size>commBuffer[1].size())
				commBuffer[1].resize(size);
			error[1]=pc.readFromNeighbor(3,commBuffer[1],numRead[1]);
			if(error[0] || error[1]) {
				printf("Error in reading from MPI neighbor %d %d.\n",error[0],error[1]);
				exit(1);
			}

			// Assign new values based on data received.
			numMoved[0]=(int)commBuffer[0][0];
			numMoved[1]=(int)commBuffer[1][0];
			numBuffer[0]=(int)commBuffer[0][1];
			numBuffer[1]=(int)commBuffer[1][1];
//			printf("Check 0 %d - %d\n",numRead[0],2+7*(numMoved[0]+numBuffer[0]));
//			printf("Check 1 %d - %d\n",numRead[1],2+7*(numMoved[1]+numBuffer[1]));
			printf("%d - Received %d %d %d %d\n",pc.getProcessNum(),numMoved[0],numBuffer[0],numMoved[1],numBuffer[1]);
//			if(pc.getProcessNum()==1) {
//				printf("0 - ");
//				for(int i=0; i<2+7*(numMoved[0]+numBuffer[0]); i++) printf("%e ",commBuffer[0][i]);
//				printf("\n");
//				printf("1 - ");
//				for(int i=0; i<2+7*(numMoved[1]+numBuffer[1]); i++) printf("%e ",commBuffer[1][i]);
//				printf("\n");
//			}
			int numBodies=pop.getNumReal();
			pop.setNumReal(pop.getNumReal()-removed.size()+numMoved[0]+numMoved[1]);
			pop.setNumBodies(pop.getNumReal()+numBuffer[0]+numBuffer[1]);

			// Fix the particles that moved.
			unsigned int copiedOver=0;		// Keeps track of how many removed
				// have been copied over by new ones coming in.
			for(int j=0; j<2; ++j) {
				for(int i=0; i<numMoved[j]; ++i) {
					if(copiedOver<removed.size()) {
						pop.setx(removed[copiedOver],commBuffer[j][2+i*7]);
						pop.sety(removed[copiedOver],commBuffer[j][3+i*7]);
						pop.setz(removed[copiedOver],commBuffer[j][4+i*7]);
						pop.setvx(removed[copiedOver],commBuffer[j][5+i*7]);
						pop.setvy(removed[copiedOver],commBuffer[j][6+i*7]);
						pop.setvz(removed[copiedOver],commBuffer[j][7+i*7]);
						pop.setRadius(removed[copiedOver],commBuffer[j][8+i*7]);
						pop.adjustAfterForce(removed[copiedOver]);
//						printf("%d - Setting real %d to ",pc.getProcessNum(),removed[copiedOver]);
//						for(int k=2; k<9; ++k) printf("%e ",commBuffer[j][k+i*7]);
//						printf("\n");
						copiedOver++;
					} else {
						pop.setx(numBodies,commBuffer[j][2+i*7]);
						pop.sety(numBodies,commBuffer[j][3+i*7]);
						pop.setz(numBodies,commBuffer[j][4+i*7]);
						pop.setvx(numBodies,commBuffer[j][5+i*7]);
						pop.setvy(numBodies,commBuffer[j][6+i*7]);
						pop.setvz(numBodies,commBuffer[j][7+i*7]);
						pop.setRadius(numBodies,commBuffer[j][8+i*7]);
						pop.adjustAfterForce(numBodies);
//						printf("%d - Setting real %d to ",pc.getProcessNum(),numBodies);
//						for(int k=2; k<9; ++k) printf("%e ",commBuffer[j][k+i*7]);
//						printf("\n");
						numBodies++;
					}
				}
			}

			// Cover extra holes if any.
			if(copiedOver<removed.size()) {
				for(unsigned int i=copiedOver; i<removed.size() && removed[i]<numBodies-1; ++i) {
					bool flag=true;
					while(flag && numBodies-1>removed[i]) {
						flag=false;
						for(unsigned int j=i+1; j<removed.size(); ++j)
							if(removed[j]==numBodies-1) {
								numBodies--;
								flag=true;
							}
					}
					if(removed[i]<numBodies-1) {
						pop.setx(removed[i],pop.getx(numBodies-1));
						pop.sety(removed[i],pop.gety(numBodies-1));
						pop.setz(removed[i],pop.getz(numBodies-1));
						pop.setvx(removed[i],pop.getvx(numBodies-1));
						pop.setvy(removed[i],pop.getvy(numBodies-1));
						pop.setvz(removed[i],pop.getvz(numBodies-1));
						pop.setX(removed[i],pop.getX(numBodies-1));
						pop.setY(removed[i],pop.getY(numBodies-1));
						pop.sete(removed[i],pop.gete(numBodies-1));
						pop.seti(removed[i],pop.geti(numBodies-1));
						pop.setPhi(removed[i],pop.getPhi(numBodies-1));
						pop.setZeta(removed[i],pop.getZeta(numBodies-1));
						pop.setRadius(removed[i],pop.getRadius(numBodies-1));
//						printf("Copying %d into %d\n",numBodies-1,removed[i]);
						numBodies--;	// Put new particles one back.
					}
				}
			}
//			printf("Critical %d==%d\n",numBodies,pop.getNumReal());

			// Add buffer particles.
			for(int j=0; j<2; ++j) {
				for(int i=numMoved[j]; i<numMoved[j]+numBuffer[j]; ++i) {
					pop.setx(numBodies,commBuffer[j][2+i*7]);
					pop.sety(numBodies,commBuffer[j][3+i*7]);
					pop.setz(numBodies,commBuffer[j][4+i*7]);
					pop.setvx(numBodies,commBuffer[j][5+i*7]);
					pop.setvy(numBodies,commBuffer[j][6+i*7]);
					pop.setvz(numBodies,commBuffer[j][7+i*7]);
					pop.setRadius(numBodies,commBuffer[j][8+i*7]);
					pop.adjustAfterForce(numBodies);
//					printf("%d - Setting buffer %d from %d to ",pc.getProcessNum(),numBodies,i);
//					for(int k=2; k<9; ++k) printf("%e ",commBuffer[j][k+i*7]);
//					printf("\n");
					numBodies++;
				}
			}

//			printf("%d - Done with BC  %d, %d\n",pc.getProcessNum(),pop.getNumReal(),pop.getNumBodies());
		}

		bool checkIfIn(double x,double y) {
			return (x>=minx && x<=maxx && y>=miny && y<=maxy);
		}

		double getArea() {
			return (maxx-minx)*(maxy-miny);
		}

		double getMinX() { return minx; }
		double getMaxX() { return maxx; }
		double getMinY() { return miny; }
		double getMaxY() { return maxy; }
		double getMinYTotal() { return minyTot; }
		double getMaxYTotal() { return maxyTot; }
	private:
		double minx,miny,maxx,maxy;

		// These values are for the whole simulation, not just on this process.
		double minxTot,maxxTot,minyTot,maxyTot;
		ProcessorCommunication &pc;
		SpacingMeasure &spaceMeasure;

		// These are the vector of particle data to send to neighbors.
		// They are encoded as numMoved,numBuffer,x1,y1,z1,vx1,vy1,vz1,x2,...
		std::vector<double> commBuffer[2];
		std::vector<int> removed;
		std::vector<double> sizeBuffer;
};


// Same as earlier, but now written to handle running in parallel.
template<class SpacingMeasure>
class ParallelSingleOrbitAzimuthal {
	public:
		ParallelSingleOrbitAzimuthal(double minX,double maxX,double minY,ProcessorCommunication &procComm,SpacingMeasure &sm):minxTot(minX),maxxTot(maxX),minyTot(minY),pc(procComm),spaceMeasure(sm) {
			minx=minxTot;
			maxx=maxxTot;
			double maxyTot=minyTot+4.0*3.141592654*GCCoords::A0*maxx;
			double size=(maxyTot-minyTot)/pc.getNumProcesses();
			miny=minyTot+pc.getProcessNum()*size;
			maxy=minyTot+(pc.getProcessNum()+1)*size;
		}

		template<class Population>
		void apply(Population &pop) {
			double bufferSize=3.0*spaceMeasure.getMinSpacing();
			miny-=2.0*GCCoords::A0*minx*pop.getTimeStep();
			maxy-=2.0*GCCoords::A0*minx*pop.getTimeStep();
			minyTot-=2.0*GCCoords::A0*minx*pop.getTimeStep();
			int numMoved[2]={0,0};
			int numBuffer[2]={0,0};
			commBuffer[0].resize(2);
			commBuffer[1].resize(2);
			removed.resize(0);
			sizeBuffer.resize(2);

			double deltay=4.0*3.141592654*GCCoords::A0*maxx;
			double tmpy=miny;
			while(tmpy<-3.141592654) tmpy+=6.283185307;
			bool wrapAround=(tmpy<-1.5*deltay || tmpy>deltay*0.5);

//			printf("%d - Applying BC %e %e %e\n",pc.getProcessNum(),miny,maxy,bufferSize);

			// Only go up to numReal.  Others are thrown away.
			for(int i=0; i<pop.getNumReal(); ++i) {
				double yTotDiff=4.0*3.141592654*GCCoords::A0*pop.getX(i);
				// Check downstream.
				if(pop.gety(i)<miny) {
					if(pc.getProcessNum()==0 && wrapAround) {
						pop.setY(i,pop.getY(i)+yTotDiff);
						pop.setCartAfterForce(i);
					}
					if(wrapAround || pc.getProcessNum()>0) {
						commBuffer[0].resize(commBuffer[0].size()+7);
						if(numBuffer[0]>0) {
							for(int j=0; j<7; ++j) {
								commBuffer[0][2+7*numMoved[0]+7*numBuffer[0]+j]=
								commBuffer[0][2+7*numMoved[0]+j];
							}
						}
						commBuffer[0][2+7*numMoved[0]]=pop.getx(i);
						commBuffer[0][3+7*numMoved[0]]=pop.gety(i);
						commBuffer[0][4+7*numMoved[0]]=pop.getz(i);
						commBuffer[0][5+7*numMoved[0]]=pop.getvx(i);
						commBuffer[0][6+7*numMoved[0]]=pop.getvy(i);
						commBuffer[0][7+7*numMoved[0]]=pop.getvz(i);
						commBuffer[0][8+7*numMoved[0]]=pop.getRadius(i);
						numMoved[0]++;
						removed.push_back(i);
					}
				} else if(pop.gety(i)<miny+bufferSize && pc.getProcessNum()!=0) {
					commBuffer[0].resize(commBuffer[0].size()+7);
					commBuffer[0][2+7*numMoved[0]+7*numBuffer[0]]=pop.getx(i);
					commBuffer[0][3+7*numMoved[0]+7*numBuffer[0]]=pop.gety(i);
					commBuffer[0][4+7*numMoved[0]+7*numBuffer[0]]=pop.getz(i);
					commBuffer[0][5+7*numMoved[0]+7*numBuffer[0]]=pop.getvx(i);
					commBuffer[0][6+7*numMoved[0]+7*numBuffer[0]]=pop.getvy(i);
					commBuffer[0][7+7*numMoved[0]+7*numBuffer[0]]=pop.getvz(i);
					commBuffer[0][8+7*numMoved[0]+7*numBuffer[0]]=pop.getRadius(i);
					numBuffer[0]++;
				} else if(pop.gety(i)>maxy || pop.gety(i)>minyTot+yTotDiff) {
				// Check upstream.
					if(pc.getProcessNum()==pc.getNumProcesses()-1 && wrapAround) {
						pop.setY(i,pop.getY(i)-yTotDiff);
						pop.setCartAfterForce(i);
					}
					if(wrapAround || pc.getProcessNum()<pc.getNumProcesses()-1) {
						commBuffer[1].resize(commBuffer[1].size()+7);
						if(numBuffer[1]>0) {
							for(int j=0; j<7; ++j) {
								commBuffer[1][2+7*numMoved[1]+7*numBuffer[1]+j]=
								commBuffer[1][2+7*numMoved[1]+j];
							}
						}
						commBuffer[1][2+7*numMoved[1]]=pop.getx(i);
						commBuffer[1][3+7*numMoved[1]]=pop.gety(i);
						commBuffer[1][4+7*numMoved[1]]=pop.getz(i);
						commBuffer[1][5+7*numMoved[1]]=pop.getvx(i);
						commBuffer[1][6+7*numMoved[1]]=pop.getvy(i);
						commBuffer[1][7+7*numMoved[1]]=pop.getvz(i);
						commBuffer[1][8+7*numMoved[1]]=pop.getRadius(i);
						numMoved[1]++;
						removed.push_back(i);
					}
				} else if(pop.gety(i)>maxy-bufferSize && pc.getProcessNum()>pc.getNumProcesses()-1) {
					commBuffer[1].resize(commBuffer[1].size()+7);
					commBuffer[1][2+7*numMoved[1]+7*numBuffer[1]]=pop.getx(i);
					commBuffer[1][3+7*numMoved[1]+7*numBuffer[1]]=pop.gety(i);
					commBuffer[1][4+7*numMoved[1]+7*numBuffer[1]]=pop.getz(i);
					commBuffer[1][5+7*numMoved[1]+7*numBuffer[1]]=pop.getvx(i);
					commBuffer[1][6+7*numMoved[1]+7*numBuffer[1]]=pop.getvy(i);
					commBuffer[1][7+7*numMoved[1]+7*numBuffer[1]]=pop.getvz(i);
					commBuffer[1][8+7*numMoved[1]+7*numBuffer[1]]=pop.getRadius(i);
					numBuffer[1]++;
				}
			}
			// Send data to neighbors.
			printf("%d - Sending in BC %d %d %d %d\n",pc.getProcessNum(),numMoved[0],numBuffer[0],numMoved[1],numBuffer[1]);
//			printf("Check 0 %d - %d\n",commBuffer[0].size(),2+7*(numMoved[0]+numBuffer[0]));
//			printf("Check 1 %d - %d\n",commBuffer[1].size(),2+7*(numMoved[1]+numBuffer[1]));
			commBuffer[0][0]=numMoved[0];
			commBuffer[0][1]=numBuffer[0];
			commBuffer[1][0]=numMoved[1];
			commBuffer[1][1]=numBuffer[1];

			// Read data from neighbors.
//			printf("%d - Reading in BC\n",pc.getProcessNum());
			int numRead[2];
			int error[2];

			// Send right, read left.
			sizeBuffer[0]=numMoved[0];
			sizeBuffer[1]=numBuffer[0];
			pc.sendToNeighbor(3,sizeBuffer);
			pc.sendToNeighbor(3,commBuffer[0]);
			error[0]=pc.readFromNeighbor(5,sizeBuffer,numRead[0]);
			unsigned int size=2+7*(int)(sizeBuffer[0]+sizeBuffer[1]);
			if(size>commBuffer[0].size())
				commBuffer[0].resize(size);
			error[0]=pc.readFromNeighbor(5,commBuffer[0],numRead[0]);

			// Send left, read right.
			sizeBuffer[0]=numMoved[1];
			sizeBuffer[1]=numBuffer[1];
			pc.sendToNeighbor(5,sizeBuffer);
			pc.sendToNeighbor(5,commBuffer[1]);
			error[1]=pc.readFromNeighbor(3,sizeBuffer,numRead[1]);
			size=2+7*(int)(sizeBuffer[0]+sizeBuffer[1]);
			if(size>commBuffer[1].size())
				commBuffer[1].resize(size);
			error[1]=pc.readFromNeighbor(3,commBuffer[1],numRead[1]);
			
			if(error[0] || error[1]) {
				printf("Error in reading from MPI neighbor %d %d.\n",error[0],error[1]);
				exit(1);
			}

			// Assign new values based on data received.
			numMoved[0]=(int)commBuffer[0][0];
			numMoved[1]=(int)commBuffer[1][0];
			numBuffer[0]=(int)commBuffer[0][1];
			numBuffer[1]=(int)commBuffer[1][1];
//			printf("Check 0 %d - %d\n",numRead[0],2+7*(numMoved[0]+numBuffer[0]));
//			printf("Check 1 %d - %d\n",numRead[1],2+7*(numMoved[1]+numBuffer[1]));
			printf("%d - Received %d %d %d %d\n",pc.getProcessNum(),numMoved[0],numBuffer[0],numMoved[1],numBuffer[1]);
			int numBodies=pop.getNumReal();
			pop.setNumReal(pop.getNumReal()-removed.size()+numMoved[0]+numMoved[1]);
			pop.setNumBodies(pop.getNumReal()+numBuffer[0]+numBuffer[1]);

			// Fix the particles that moved.
			unsigned int copiedOver=0;		// Keeps track of how many removed
				// have been copied over by new ones coming in.
			for(int j=0; j<2; ++j) {
				for(int i=0; i<numMoved[j]; ++i) {
					if(copiedOver<removed.size()) {
						pop.setx(removed[copiedOver],commBuffer[j][2+i*7]);
						pop.sety(removed[copiedOver],commBuffer[j][3+i*7]);
						pop.setz(removed[copiedOver],commBuffer[j][4+i*7]);
						pop.setvx(removed[copiedOver],commBuffer[j][5+i*7]);
						pop.setvy(removed[copiedOver],commBuffer[j][6+i*7]);
						pop.setvz(removed[copiedOver],commBuffer[j][7+i*7]);
						pop.setRadius(removed[copiedOver],commBuffer[j][8+i*7]);
						pop.adjustAfterForce(removed[copiedOver]);
//						printf("%d - Setting real %d to ",pc.getProcessNum(),removed[copiedOver]);
//						for(int k=2; k<9; ++k) printf("%e ",commBuffer[j][k+i*7]);
//						printf("\n");
						copiedOver++;
					} else {
						pop.setx(numBodies,commBuffer[j][2+i*7]);
						pop.sety(numBodies,commBuffer[j][3+i*7]);
						pop.setz(numBodies,commBuffer[j][4+i*7]);
						pop.setvx(numBodies,commBuffer[j][5+i*7]);
						pop.setvy(numBodies,commBuffer[j][6+i*7]);
						pop.setvz(numBodies,commBuffer[j][7+i*7]);
						pop.setRadius(numBodies,commBuffer[j][8+i*7]);
						pop.adjustAfterForce(numBodies);
//						printf("%d - Setting real %d to ",pc.getProcessNum(),numBodies);
//						for(int k=2; k<9; ++k) printf("%e ",commBuffer[j][k+i*7]);
//						printf("\n");
						numBodies++;
					}
				}
			}

			// Cover extra holes if any.
			if(copiedOver<removed.size()) {
				for(unsigned int i=copiedOver; i<removed.size() && removed[i]<numBodies-1; ++i) {
					bool flag=true;
					while(flag && numBodies-1>removed[i]) {
						flag=false;
						for(unsigned int j=i+1; j<removed.size(); ++j)
							if(removed[j]==numBodies-1) {
								numBodies--;
								flag=true;
							}
					}
					if(removed[i]<numBodies-1) {
						pop.setx(removed[i],pop.getx(numBodies-1));
						pop.sety(removed[i],pop.gety(numBodies-1));
						pop.setz(removed[i],pop.getz(numBodies-1));
						pop.setvx(removed[i],pop.getvx(numBodies-1));
						pop.setvy(removed[i],pop.getvy(numBodies-1));
						pop.setvz(removed[i],pop.getvz(numBodies-1));
						pop.setX(removed[i],pop.getX(numBodies-1));
						pop.setY(removed[i],pop.getY(numBodies-1));
						pop.sete(removed[i],pop.gete(numBodies-1));
						pop.seti(removed[i],pop.geti(numBodies-1));
						pop.setPhi(removed[i],pop.getPhi(numBodies-1));
						pop.setZeta(removed[i],pop.getZeta(numBodies-1));
						pop.setRadius(removed[i],pop.getRadius(numBodies-1));
//						printf("Copying %d into %d\n",numBodies-1,removed[i]);
						numBodies--;	// Put new particles one back.
					}
				}
			}
//			printf("Critical %d==%d\n",numBodies,pop.getNumReal());

			// Add buffer particles.
			for(int j=0; j<2; ++j) {
				for(int i=numMoved[j]; i<numMoved[j]+numBuffer[j]; ++i) {
					pop.setx(numBodies,commBuffer[j][2+i*7]);
					pop.sety(numBodies,commBuffer[j][3+i*7]);
					pop.setz(numBodies,commBuffer[j][4+i*7]);
					pop.setvx(numBodies,commBuffer[j][5+i*7]);
					pop.setvy(numBodies,commBuffer[j][6+i*7]);
					pop.setvz(numBodies,commBuffer[j][7+i*7]);
					pop.setRadius(numBodies,commBuffer[j][8+i*7]);
					pop.adjustAfterForce(numBodies);
//					printf("%d - Setting buffer %d from %d to ",pc.getProcessNum(),numBodies,i);
//					for(int k=2; k<9; ++k) printf("%e ",commBuffer[j][k+i*7]);
//					printf("\n");
					numBodies++;
				}
			}

//			printf("%d - Done with BC  %d, %d\n",pc.getProcessNum(),pop.getNumReal(),pop.getNumBodies());
		}

		bool checkIfIn(double x,double y) {
			return (x>=minx && x<=maxx && y>=miny && y<=maxy && y<=minyTot+4.0*3.141592654*GCCoords::A0*x);
		}

		double getArea() {
			return (getMaxX()-getMinX())*(getMaxY()-getMinY());
/*			double basey=minyTot+4.0*3.141592654*GCCoords::A0*minx;
			if(maxy<=basey) {
				return (maxx-minx)*(maxy-miny);
			} else {
				double ret=0.0;
				double basex,topx;
				if(miny<basey) {
					ret+=(maxx-minx)*(basey-miny);
					basex=minx;
				} else {
					basex=(miny-minyTot)/(4.0*3.141592654*GCCoords::A0);
				}
				double smally=(basey>miny)?basey:miny;
				topx=(maxy-minyTot)/(4.0*3.141592654*GCCoords::A0);
				return ret+(maxy-smally)*0.5*((maxx-basex)+(maxx-topx));
			}
*/		}

		double getMinX() { return minx; }
		double getMaxX() { return maxx; }
		double getMinY() { return miny; }
		double getMaxY() { return maxy; }
		double getMinYTotal() { return minyTot; }
		double getMaxYTotal() { return minyTot+4.0*3.141592654*GCCoords::A0*maxx; }
	private:
		double minx,miny,maxx,maxy;

		// These values are for the whole simulation, not just on this process.
		double minxTot,maxxTot,minyTot;
		ProcessorCommunication &pc;
		SpacingMeasure &spaceMeasure;

		// These are the vector of particle data to send to neighbors.
		// They are encoded as numMoved,numBuffer,x1,y1,z1,vx1,vy1,vz1,x2,...
		std::vector<double> commBuffer[2];
		std::vector<int> removed;
		std::vector<double> sizeBuffer;
};

template<class SpacingMeasure>
class ParallelFixedPeriodic {
	public:
		ParallelFixedPeriodic(double minX,double maxX,double minY,double maxY,ProcessorCommunication &procComm,SpacingMeasure &sm,bool zero=false):minxTot(minX),maxxTot(maxX),minyTot(minY),maxyTot(maxY),pc(procComm),spaceMeasure(sm),zeroeOnWrap(zero) {
			minx=minxTot;
			maxx=maxxTot;
			double size=(maxyTot-minyTot)/pc.getNumProcesses();
			miny=minyTot+pc.getProcessNum()*size;
			maxy=minyTot+(pc.getProcessNum()+1)*size;
		}

		template<class Population>
		void apply(Population &pop) {
			double bufferSize=3.0*spaceMeasure.getMinSpacing();
			int numMoved[2]={0,0};
			int numBuffer[2]={0,0};
			commBuffer[0].resize(2);
			commBuffer[1].resize(2);
			removed.resize(0);
			sizeBuffer.resize(2);

//			printf("%d - Applying BC %e %e %e %d\n",pc.getProcessNum(),miny,maxy,bufferSize,pop.getNumReal());

			// Only go up to numReal.  Others are thrown away.
			for(int i=0; i<pop.getNumReal(); ++i) {
				// Check downstream.
//				if(i%10000==0) printf("i[%d]=%e\n",i,pop.gety(i));
				if(pop.gety(i)<miny) {
					if(pc.getProcessNum()==0) {
						pop.setY(i,pop.getY(i)+(maxyTot-minyTot));
						if(zeroeOnWrap) pop.sete(i,0);
						pop.setCartAfterForce(i);
					}
					commBuffer[0].resize(commBuffer[0].size()+7);
					if(numBuffer[0]>0) {
						for(int j=0; j<7; ++j) {
							commBuffer[0][2+7*numMoved[0]+7*numBuffer[0]+j]=
							commBuffer[0][2+7*numMoved[0]+j];
						}
					}
					commBuffer[0][2+7*numMoved[0]]=pop.getx(i);
					commBuffer[0][3+7*numMoved[0]]=pop.gety(i);
					commBuffer[0][4+7*numMoved[0]]=pop.getz(i);
					commBuffer[0][5+7*numMoved[0]]=pop.getvx(i);
					commBuffer[0][6+7*numMoved[0]]=pop.getvy(i);
					commBuffer[0][7+7*numMoved[0]]=pop.getvz(i);
					commBuffer[0][8+7*numMoved[0]]=pop.getRadius(i);
					numMoved[0]++;
					removed.push_back(i);
				} else if(pop.gety(i)<miny+bufferSize && pc.getProcessNum()>0) {
					commBuffer[0].resize(commBuffer[0].size()+7);
					commBuffer[0][2+7*numMoved[0]+7*numBuffer[0]]=pop.getx(i);
					commBuffer[0][3+7*numMoved[0]+7*numBuffer[0]]=pop.gety(i);
					commBuffer[0][4+7*numMoved[0]+7*numBuffer[0]]=pop.getz(i);
					commBuffer[0][5+7*numMoved[0]+7*numBuffer[0]]=pop.getvx(i);
					commBuffer[0][6+7*numMoved[0]+7*numBuffer[0]]=pop.getvy(i);
					commBuffer[0][7+7*numMoved[0]+7*numBuffer[0]]=pop.getvz(i);
					commBuffer[0][8+7*numMoved[0]+7*numBuffer[0]]=pop.getRadius(i);
					numBuffer[0]++;
				} else if(pop.gety(i)>maxy) {
				// Check upstream.
					if(pc.getProcessNum()==pc.getNumProcesses()-1) {
						pop.setY(i,pop.getY(i)-(maxyTot-minyTot));
						pop.setCartAfterForce(i);
					}
					commBuffer[1].resize(commBuffer[1].size()+7);
					if(numBuffer[1]>0) {
						for(int j=0; j<7; ++j) {
							commBuffer[1][2+7*numMoved[1]+7*numBuffer[1]+j]=
							commBuffer[1][2+7*numMoved[1]+j];
						}
					}
					commBuffer[1][2+7*numMoved[1]]=pop.getx(i);
					commBuffer[1][3+7*numMoved[1]]=pop.gety(i);
					commBuffer[1][4+7*numMoved[1]]=pop.getz(i);
					commBuffer[1][5+7*numMoved[1]]=pop.getvx(i);
					commBuffer[1][6+7*numMoved[1]]=pop.getvy(i);
					commBuffer[1][7+7*numMoved[1]]=pop.getvz(i);
					commBuffer[1][8+7*numMoved[1]]=pop.getRadius(i);
					numMoved[1]++;
					removed.push_back(i);
				} else if(pop.gety(i)>maxy-bufferSize && pc.getProcessNum()<pc.getNumProcesses()-1) {
					commBuffer[1].resize(commBuffer[1].size()+7);
					commBuffer[1][2+7*numMoved[1]+7*numBuffer[1]]=pop.getx(i);
					commBuffer[1][3+7*numMoved[1]+7*numBuffer[1]]=pop.gety(i);
					commBuffer[1][4+7*numMoved[1]+7*numBuffer[1]]=pop.getz(i);
					commBuffer[1][5+7*numMoved[1]+7*numBuffer[1]]=pop.getvx(i);
					commBuffer[1][6+7*numMoved[1]+7*numBuffer[1]]=pop.getvy(i);
					commBuffer[1][7+7*numMoved[1]+7*numBuffer[1]]=pop.getvz(i);
					commBuffer[1][8+7*numMoved[1]+7*numBuffer[1]]=pop.getRadius(i);
					numBuffer[1]++;
				}

				// Check inward if neighbors that direction.
				// !!! Not written yet.

				// Check outward if neighbors that direction.
				// !!! Not written yet.
			}
			// Send data to neighbors.
			printf("%d - Sending in BC %d %d %d %d\n",pc.getProcessNum(),numMoved[0],numBuffer[0],numMoved[1],numBuffer[1]);
//			printf("Check 0 %d - %d\n",commBuffer[0].size(),2+7*(numMoved[0]+numBuffer[0]));
//			printf("Check 1 %d - %d\n",commBuffer[1].size(),2+7*(numMoved[1]+numBuffer[1]));
//			if(pc.getProcessNum()==1) {
//				printf("Sending 0 - ");
//				for(int i=0; i<commBuffer[0].size(); i++) printf("%e ",commBuffer[0][i]);
//				printf("\n");
//				printf("Sending 1 - ");
//				for(int i=0; i<2+7*(numMoved[1]+numBuffer[1]); i++) printf("%e ",commBuffer[1][i]);
//				printf("\n");
//			}
			commBuffer[0][0]=numMoved[0];
			commBuffer[0][1]=numBuffer[0];
			commBuffer[1][0]=numMoved[1];
			commBuffer[1][1]=numBuffer[1];

			// Read data from neighbors.
//			printf("%d - Reading in BC\n",pc.getProcessNum());
			int numRead[2];
			int error[2];

			// Send right, read left.
			if(pc.getProcessNum()==0) {
				printf("%d - Send right\n",pc.getProcessNum());
				sizeBuffer[0]=numMoved[0];
				sizeBuffer[1]=numBuffer[0];
				pc.sendToNeighbor(3,sizeBuffer);
				pc.sendToNeighbor(3,commBuffer[0]);
				printf("%d - Read left\n",pc.getProcessNum());
				error[0]=pc.readFromNeighbor(5,sizeBuffer,numRead[0]);
				unsigned int size=2+7*(int)(sizeBuffer[0]+sizeBuffer[1]);
				if(size>commBuffer[0].size())
					commBuffer[0].resize(size);
				error[0]=pc.readFromNeighbor(5,commBuffer[0],numRead[0]);
			} else {
				printf("%d - Read left\n",pc.getProcessNum());
				printf("%d - Read size\n",pc.getProcessNum());
				error[0]=pc.readFromNeighbor(5,sizeBuffer,numRead[0]);
				unsigned int size=2+7*(int)(sizeBuffer[0]+sizeBuffer[1]);
				printf("Resize to %d\n",size);
				if(size>recBuffer.size())
					recBuffer.resize(size);
				printf("%d - Read bodies\n",pc.getProcessNum());
				error[0]=pc.readFromNeighbor(5,recBuffer,numRead[0]);
				printf("%d - Send right\n",pc.getProcessNum());
				sizeBuffer[0]=numMoved[0];
				sizeBuffer[1]=numBuffer[0];
				pc.sendToNeighbor(3,sizeBuffer);
				pc.sendToNeighbor(3,commBuffer[0]);
				commBuffer[0]=recBuffer;
			}

			// Send left, read right.
			if(pc.getProcessNum()==0) {
				printf("%d - Send left\n",pc.getProcessNum());
				sizeBuffer[0]=numMoved[1];
				sizeBuffer[1]=numBuffer[1];
				pc.sendToNeighbor(5,sizeBuffer);
				pc.sendToNeighbor(5,commBuffer[1]);
				printf("%d - Read right\n",pc.getProcessNum());
				error[1]=pc.readFromNeighbor(3,sizeBuffer,numRead[1]);
				unsigned int size=2+7*(int)(sizeBuffer[0]+sizeBuffer[1]);
				if(size>commBuffer[1].size())
					commBuffer[1].resize(size);
				error[1]=pc.readFromNeighbor(3,commBuffer[1],numRead[1]);
			} else {
				printf("%d - Read right\n",pc.getProcessNum());
				error[1]=pc.readFromNeighbor(3,sizeBuffer,numRead[1]);
				unsigned int size=2+7*(int)(sizeBuffer[0]+sizeBuffer[1]);
				if(size>recBuffer.size())
					recBuffer.resize(size);
				error[1]=pc.readFromNeighbor(3,recBuffer,numRead[1]);
				printf("%d - Send left\n",pc.getProcessNum());
				sizeBuffer[0]=numMoved[1];
				sizeBuffer[1]=numBuffer[1];
				pc.sendToNeighbor(5,sizeBuffer);
				pc.sendToNeighbor(5,commBuffer[1]);
				commBuffer[1]=recBuffer;
			}

			if(error[0] || error[1]) {
				printf("Error in reading from MPI neighbor %d %d.\n",error[0],error[1]);
				exit(1);
			}

			// Assign new values based on data received.
			numMoved[0]=(int)commBuffer[0][0];
			numMoved[1]=(int)commBuffer[1][0];
			numBuffer[0]=(int)commBuffer[0][1];
			numBuffer[1]=(int)commBuffer[1][1];
//			printf("Check 0 %d - %d\n",numRead[0],2+7*(numMoved[0]+numBuffer[0]));
//			printf("Check 1 %d - %d\n",numRead[1],2+7*(numMoved[1]+numBuffer[1]));
			printf("%d - Received %d %d %d %d\n",pc.getProcessNum(),numMoved[0],numBuffer[0],numMoved[1],numBuffer[1]);
//			if(pc.getProcessNum()==1) {
//				printf("0 - ");
//				for(int i=0; i<2+7*(numMoved[0]+numBuffer[0]); i++) printf("%e ",commBuffer[0][i]);
//				printf("\n");
//				printf("1 - ");
//				for(int i=0; i<2+7*(numMoved[1]+numBuffer[1]); i++) printf("%e ",commBuffer[1][i]);
//				printf("\n");
//			}
			int numBodies=pop.getNumReal();
			pop.setNumReal(pop.getNumReal()-removed.size()+numMoved[0]+numMoved[1]);
			pop.setNumBodies(pop.getNumReal()+numBuffer[0]+numBuffer[1]);

			// Fix the particles that moved.
			unsigned int copiedOver=0;		// Keeps track of how many removed
				// have been copied over by new ones coming in.
			for(int j=0; j<2; ++j) {
				for(int i=0; i<numMoved[j]; ++i) {
					if(copiedOver<removed.size()) {
						pop.setx(removed[copiedOver],commBuffer[j][2+i*7]);
						pop.sety(removed[copiedOver],commBuffer[j][3+i*7]);
						pop.setz(removed[copiedOver],commBuffer[j][4+i*7]);
						pop.setvx(removed[copiedOver],commBuffer[j][5+i*7]);
						pop.setvy(removed[copiedOver],commBuffer[j][6+i*7]);
						pop.setvz(removed[copiedOver],commBuffer[j][7+i*7]);
						pop.setRadius(removed[copiedOver],commBuffer[j][8+i*7]);
						pop.adjustAfterForce(removed[copiedOver]);
//						printf("%d - Setting real %d to ",pc.getProcessNum(),removed[copiedOver]);
//						for(int k=2; k<9; ++k) printf("%e ",commBuffer[j][k+i*7]);
//						printf("\n");
						copiedOver++;
					} else {
						pop.setx(numBodies,commBuffer[j][2+i*7]);
						pop.sety(numBodies,commBuffer[j][3+i*7]);
						pop.setz(numBodies,commBuffer[j][4+i*7]);
						pop.setvx(numBodies,commBuffer[j][5+i*7]);
						pop.setvy(numBodies,commBuffer[j][6+i*7]);
						pop.setvz(numBodies,commBuffer[j][7+i*7]);
						pop.setRadius(numBodies,commBuffer[j][8+i*7]);
						pop.adjustAfterForce(numBodies);
//						printf("%d - Setting real %d to ",pc.getProcessNum(),numBodies);
//						for(int k=2; k<9; ++k) printf("%e ",commBuffer[j][k+i*7]);
//						printf("\n");
						numBodies++;
					}
				}
			}

			// Cover extra holes if any.
			if(copiedOver<removed.size()) {
				for(unsigned int i=copiedOver; i<removed.size() && removed[i]<numBodies-1; ++i) {
					bool flag=true;
					while(flag && numBodies-1>removed[i]) {
						flag=false;
						for(unsigned int j=i+1; j<removed.size(); ++j)
							if(removed[j]==numBodies-1) {
								numBodies--;
								flag=true;
							}
					}
					if(removed[i]<numBodies-1) {
						pop.setx(removed[i],pop.getx(numBodies-1));
						pop.sety(removed[i],pop.gety(numBodies-1));
						pop.setz(removed[i],pop.getz(numBodies-1));
						pop.setvx(removed[i],pop.getvx(numBodies-1));
						pop.setvy(removed[i],pop.getvy(numBodies-1));
						pop.setvz(removed[i],pop.getvz(numBodies-1));
						pop.setX(removed[i],pop.getX(numBodies-1));
						pop.setY(removed[i],pop.getY(numBodies-1));
						pop.sete(removed[i],pop.gete(numBodies-1));
						pop.seti(removed[i],pop.geti(numBodies-1));
						pop.setPhi(removed[i],pop.getPhi(numBodies-1));
						pop.setZeta(removed[i],pop.getZeta(numBodies-1));
						pop.setRadius(removed[i],pop.getRadius(numBodies-1));
//						printf("Copying %d into %d\n",numBodies-1,removed[i]);
						numBodies--;	// Put new particles one back.
					}
				}
			}
			printf("Critical %d==%d\n",numBodies,pop.getNumReal());

			// Add buffer particles.
			for(int j=0; j<2; ++j) {
				for(int i=numMoved[j]; i<numMoved[j]+numBuffer[j]; ++i) {
					pop.setx(numBodies,commBuffer[j][2+i*7]);
					pop.sety(numBodies,commBuffer[j][3+i*7]);
					pop.setz(numBodies,commBuffer[j][4+i*7]);
					pop.setvx(numBodies,commBuffer[j][5+i*7]);
					pop.setvy(numBodies,commBuffer[j][6+i*7]);
					pop.setvz(numBodies,commBuffer[j][7+i*7]);
					pop.setRadius(numBodies,commBuffer[j][8+i*7]);
					pop.adjustAfterForce(numBodies);
//					printf("%d - Setting buffer %d from %d to ",pc.getProcessNum(),numBodies,i);
//					for(int k=2; k<9; ++k) printf("%e ",commBuffer[j][k+i*7]);
//					printf("\n");
					numBodies++;
				}
			}

			printf("%d - Done with BC  %d, %d\n",pc.getProcessNum(),pop.getNumReal(),pop.getNumBodies());
		}

		bool checkIfIn(double x,double y) {
			return (x>=minx && x<=maxx && y>=miny && y<=maxy);
		}

		double getArea() {
			return (maxx-minx)*(maxy-miny);
		}

		double getMinX() { return minx; }
		double getMaxX() { return maxx; }
		double getMinY() { return miny; }
		double getMaxY() { return maxy; }
		double getMinYTotal() { return minyTot; }
		double getMaxYTotal() { return maxyTot; }
	private:
		double minx,miny,maxx,maxy;

		// These values are for the whole simulation, not just on this process.
		double minxTot,maxxTot,minyTot,maxyTot;
		ProcessorCommunication &pc;
		SpacingMeasure &spaceMeasure;
		bool zeroeOnWrap;

		// These are the vector of particle data to send to neighbors.
		// They are encoded as numMoved,numBuffer,x1,y1,z1,vx1,vy1,vz1,x2,...
		std::vector<double> commBuffer[2];
		std::vector<double> recBuffer;
		std::vector<int> removed;
		std::vector<double> sizeBuffer;
};

template<class SpacingMeasure>
class ParallelShearedPeriodic {
	public:
		ParallelShearedPeriodic(double sizeX,double sizeY,double sizeZ,double sh,ProcessorCommunication &procComm,SpacingMeasure &sm):minxTot(-sizeX),maxxTot(sizeX),minyTot(-sizeY),maxyTot(sizeY),sizez(sizeZ),shear(sh),pc(procComm),spaceMeasure(sm) {
			minx=minxTot;
			maxx=maxxTot;
			double size=(maxyTot-minyTot)/pc.getNumProcesses();
			miny=minyTot+pc.getProcessNum()*size;
			maxy=minyTot+(pc.getProcessNum()+1)*size;
		}

		template<class Population>
		void apply(Population &pop) {
			double bufferSize=3.0*spaceMeasure.getMinSpacing();
			int numMoved[2]={0,0};
			int numBuffer[2]={0,0};
			commBuffer[0].resize(2);
			commBuffer[1].resize(2);
			removed.resize(0);
			sizeBuffer.resize(2);

//			printf("%d - Applying BC %e %e %e %d\n",pc.getProcessNum(),miny,maxy,bufferSize,pop.getNumReal());

			// Only go up to numReal.  Others are thrown away.
			for(int i=0; i<pop.getNumReal(); ++i) {
				while(pop.getx(i)<minx) {
					pop.setx(i,pop.getx(i)+(maxx-minx));
				} 
				while(pop.getx(i)>maxx) {
					pop.setx(i,pop.getx(i)-(maxx-minx));
				} 
				while(pop.getz(i)<-sizez) {
					pop.setz(i,pop.getz(i)+(2.0*sizez));
				} 
				while(pop.getz(i)>sizez) {
					pop.setz(i,pop.getz(i)-(2.0*sizez));
				} 
				// Check downstream.
//				if(i%10000==0) printf("i[%d]=%e\n",i,pop.gety(i));
				if(pop.gety(i)<miny) {
					if(pc.getProcessNum()==0) {
						pop.sety(i,pop.gety(i)+(maxyTot-minyTot));
						pop.setvx(i,pop.getvx(i)+shear);
					}
					commBuffer[0].resize(commBuffer[0].size()+7);
					if(numBuffer[0]>0) {
						for(int j=0; j<7; ++j) {
							commBuffer[0][2+7*numMoved[0]+7*numBuffer[0]+j]=
							commBuffer[0][2+7*numMoved[0]+j];
						}
					}
					commBuffer[0][2+7*numMoved[0]]=pop.getx(i);
					commBuffer[0][3+7*numMoved[0]]=pop.gety(i);
					commBuffer[0][4+7*numMoved[0]]=pop.getz(i);
					commBuffer[0][5+7*numMoved[0]]=pop.getvx(i);
					commBuffer[0][6+7*numMoved[0]]=pop.getvy(i);
					commBuffer[0][7+7*numMoved[0]]=pop.getvz(i);
					commBuffer[0][8+7*numMoved[0]]=pop.getRadius(i);
					numMoved[0]++;
					removed.push_back(i);
				} else if(pop.gety(i)<miny+bufferSize && pc.getProcessNum()>0) {
					commBuffer[0].resize(commBuffer[0].size()+7);
					commBuffer[0][2+7*numMoved[0]+7*numBuffer[0]]=pop.getx(i);
					commBuffer[0][3+7*numMoved[0]+7*numBuffer[0]]=pop.gety(i);
					commBuffer[0][4+7*numMoved[0]+7*numBuffer[0]]=pop.getz(i);
					commBuffer[0][5+7*numMoved[0]+7*numBuffer[0]]=pop.getvx(i);
					commBuffer[0][6+7*numMoved[0]+7*numBuffer[0]]=pop.getvy(i);
					commBuffer[0][7+7*numMoved[0]+7*numBuffer[0]]=pop.getvz(i);
					commBuffer[0][8+7*numMoved[0]+7*numBuffer[0]]=pop.getRadius(i);
					numBuffer[0]++;
				} else if(pop.gety(i)>maxy) {
				// Check upstream.
					if(pc.getProcessNum()==pc.getNumProcesses()-1) {
						pop.sety(i,pop.gety(i)-(maxyTot-minyTot));
						pop.setvx(i,pop.getvx(i)-shear);
					}
					commBuffer[1].resize(commBuffer[1].size()+7);
					if(numBuffer[1]>0) {
						for(int j=0; j<7; ++j) {
							commBuffer[1][2+7*numMoved[1]+7*numBuffer[1]+j]=
							commBuffer[1][2+7*numMoved[1]+j];
						}
					}
					commBuffer[1][2+7*numMoved[1]]=pop.getx(i);
					commBuffer[1][3+7*numMoved[1]]=pop.gety(i);
					commBuffer[1][4+7*numMoved[1]]=pop.getz(i);
					commBuffer[1][5+7*numMoved[1]]=pop.getvx(i);
					commBuffer[1][6+7*numMoved[1]]=pop.getvy(i);
					commBuffer[1][7+7*numMoved[1]]=pop.getvz(i);
					commBuffer[1][8+7*numMoved[1]]=pop.getRadius(i);
					numMoved[1]++;
					removed.push_back(i);
				} else if(pop.gety(i)>maxy-bufferSize && pc.getProcessNum()<pc.getNumProcesses()-1) {
					commBuffer[1].resize(commBuffer[1].size()+7);
					commBuffer[1][2+7*numMoved[1]+7*numBuffer[1]]=pop.getx(i);
					commBuffer[1][3+7*numMoved[1]+7*numBuffer[1]]=pop.gety(i);
					commBuffer[1][4+7*numMoved[1]+7*numBuffer[1]]=pop.getz(i);
					commBuffer[1][5+7*numMoved[1]+7*numBuffer[1]]=pop.getvx(i);
					commBuffer[1][6+7*numMoved[1]+7*numBuffer[1]]=pop.getvy(i);
					commBuffer[1][7+7*numMoved[1]+7*numBuffer[1]]=pop.getvz(i);
					commBuffer[1][8+7*numMoved[1]+7*numBuffer[1]]=pop.getRadius(i);
					numBuffer[1]++;
				}

				// Check inward if neighbors that direction.
				// !!! Not written yet.

				// Check outward if neighbors that direction.
				// !!! Not written yet.
			}
			// Send data to neighbors.
			printf("%d - Sending in BC %d %d %d %d\n",pc.getProcessNum(),numMoved[0],numBuffer[0],numMoved[1],numBuffer[1]);
//			printf("Check 0 %d - %d\n",commBuffer[0].size(),2+7*(numMoved[0]+numBuffer[0]));
//			printf("Check 1 %d - %d\n",commBuffer[1].size(),2+7*(numMoved[1]+numBuffer[1]));
//			if(pc.getProcessNum()==1) {
//				printf("Sending 0 - ");
//				for(int i=0; i<commBuffer[0].size(); i++) printf("%e ",commBuffer[0][i]);
//				printf("\n");
//				printf("Sending 1 - ");
//				for(int i=0; i<2+7*(numMoved[1]+numBuffer[1]); i++) printf("%e ",commBuffer[1][i]);
//				printf("\n");
//			}
			commBuffer[0][0]=numMoved[0];
			commBuffer[0][1]=numBuffer[0];
			commBuffer[1][0]=numMoved[1];
			commBuffer[1][1]=numBuffer[1];

			// Read data from neighbors.
//			printf("%d - Reading in BC\n",pc.getProcessNum());
			int numRead[2];
			int error[2];

			// Send right, read left.
			sizeBuffer[0]=numMoved[0];
			sizeBuffer[1]=numBuffer[0];
			pc.sendToNeighbor(3,sizeBuffer);
			pc.sendToNeighbor(3,commBuffer[0]);
			error[0]=pc.readFromNeighbor(5,sizeBuffer,numRead[0]);
			unsigned int size=2+7*(int)(sizeBuffer[0]+sizeBuffer[1]);
			if(size>commBuffer[0].size())
				commBuffer[0].resize(size);
			error[0]=pc.readFromNeighbor(5,commBuffer[0],numRead[0]);

			// Send left, read right.
			sizeBuffer[0]=numMoved[1];
			sizeBuffer[1]=numBuffer[1];
			pc.sendToNeighbor(5,sizeBuffer);
			pc.sendToNeighbor(5,commBuffer[1]);
			error[1]=pc.readFromNeighbor(3,sizeBuffer,numRead[1]);
			size=2+7*(int)(sizeBuffer[0]+sizeBuffer[1]);
			if(size>commBuffer[1].size())
				commBuffer[1].resize(size);
			error[1]=pc.readFromNeighbor(3,commBuffer[1],numRead[1]);
			if(error[0] || error[1]) {
				printf("Error in reading from MPI neighbor %d %d.\n",error[0],error[1]);
				exit(1);
			}

			// Assign new values based on data received.
			numMoved[0]=(int)commBuffer[0][0];
			numMoved[1]=(int)commBuffer[1][0];
			numBuffer[0]=(int)commBuffer[0][1];
			numBuffer[1]=(int)commBuffer[1][1];
//			printf("Check 0 %d - %d\n",numRead[0],2+7*(numMoved[0]+numBuffer[0]));
//			printf("Check 1 %d - %d\n",numRead[1],2+7*(numMoved[1]+numBuffer[1]));
//			printf("%d - Received %d %d %d %d\n",pc.getProcessNum(),numMoved[0],numBuffer[0],numMoved[1],numBuffer[1]);
//			if(pc.getProcessNum()==1) {
//				printf("0 - ");
//				for(int i=0; i<2+7*(numMoved[0]+numBuffer[0]); i++) printf("%e ",commBuffer[0][i]);
//				printf("\n");
//				printf("1 - ");
//				for(int i=0; i<2+7*(numMoved[1]+numBuffer[1]); i++) printf("%e ",commBuffer[1][i]);
//				printf("\n");
//			}
			int numBodies=pop.getNumReal();
			pop.setNumReal(pop.getNumReal()-removed.size()+numMoved[0]+numMoved[1]);
			pop.setNumBodies(pop.getNumReal()+numBuffer[0]+numBuffer[1]);

			// Fix the particles that moved.
			unsigned int copiedOver=0;		// Keeps track of how many removed
				// have been copied over by new ones coming in.
			for(int j=0; j<2; ++j) {
				for(int i=0; i<numMoved[j]; ++i) {
					if(copiedOver<removed.size()) {
						pop.setx(removed[copiedOver],commBuffer[j][2+i*7]);
						pop.sety(removed[copiedOver],commBuffer[j][3+i*7]);
						pop.setz(removed[copiedOver],commBuffer[j][4+i*7]);
						pop.setvx(removed[copiedOver],commBuffer[j][5+i*7]);
						pop.setvy(removed[copiedOver],commBuffer[j][6+i*7]);
						pop.setvz(removed[copiedOver],commBuffer[j][7+i*7]);
						pop.setRadius(removed[copiedOver],commBuffer[j][8+i*7]);
//						printf("%d - Setting real %d to ",pc.getProcessNum(),removed[copiedOver]);
//						for(int k=2; k<9; ++k) printf("%e ",commBuffer[j][k+i*7]);
//						printf("\n");
						copiedOver++;
					} else {
						pop.setx(numBodies,commBuffer[j][2+i*7]);
						pop.sety(numBodies,commBuffer[j][3+i*7]);
						pop.setz(numBodies,commBuffer[j][4+i*7]);
						pop.setvx(numBodies,commBuffer[j][5+i*7]);
						pop.setvy(numBodies,commBuffer[j][6+i*7]);
						pop.setvz(numBodies,commBuffer[j][7+i*7]);
						pop.setRadius(numBodies,commBuffer[j][8+i*7]);
//						printf("%d - Setting real %d to ",pc.getProcessNum(),numBodies);
//						for(int k=2; k<9; ++k) printf("%e ",commBuffer[j][k+i*7]);
//						printf("\n");
						numBodies++;
					}
				}
			}

			// Cover extra holes if any.
			if(copiedOver<removed.size()) {
				for(unsigned int i=copiedOver; i<removed.size() && removed[i]<numBodies-1; ++i) {
					bool flag=true;
					while(flag && numBodies-1>removed[i]) {
						flag=false;
						for(unsigned int j=i+1; j<removed.size(); ++j)
							if(removed[j]==numBodies-1) {
								numBodies--;
								flag=true;
							}
					}
					if(removed[i]<numBodies-1) {
						pop.setx(removed[i],pop.getx(numBodies-1));
						pop.sety(removed[i],pop.gety(numBodies-1));
						pop.setz(removed[i],pop.getz(numBodies-1));
						pop.setvx(removed[i],pop.getvx(numBodies-1));
						pop.setvy(removed[i],pop.getvy(numBodies-1));
						pop.setvz(removed[i],pop.getvz(numBodies-1));
						pop.setRadius(removed[i],pop.getRadius(numBodies-1));
//						printf("Copying %d into %d\n",numBodies-1,removed[i]);
						numBodies--;	// Put new particles one back.
					}
				}
			}
//			printf("Critical %d==%d\n",numBodies,pop.getNumReal());

			// Add buffer particles.
			for(int j=0; j<2; ++j) {
				for(int i=numMoved[j]; i<numMoved[j]+numBuffer[j]; ++i) {
					pop.setx(numBodies,commBuffer[j][2+i*7]);
					pop.sety(numBodies,commBuffer[j][3+i*7]);
					pop.setz(numBodies,commBuffer[j][4+i*7]);
					pop.setvx(numBodies,commBuffer[j][5+i*7]);
					pop.setvy(numBodies,commBuffer[j][6+i*7]);
					pop.setvz(numBodies,commBuffer[j][7+i*7]);
					pop.setRadius(numBodies,commBuffer[j][8+i*7]);
//					printf("%d - Setting buffer %d from %d to ",pc.getProcessNum(),numBodies,i);
//					for(int k=2; k<9; ++k) printf("%e ",commBuffer[j][k+i*7]);
//					printf("\n");
					numBodies++;
				}
			}

//			printf("%d - Done with BC  %d, %d\n",pc.getProcessNum(),pop.getNumReal(),pop.getNumBodies());
		}

		bool checkIfIn(double x,double y) {
			return (x>=minx && x<=maxx && y>=miny && y<=maxy);
		}

		double getArea() {
			return (maxx-minx)*(maxy-miny);
		}

		double getMinX() { return minx; }
		double getMaxX() { return maxx; }
		double getMinY() { return miny; }
		double getMaxY() { return maxy; }
		double getMinZ() { return -sizez; }
		double getMaxZ() { return sizez; }
		double getMinYTotal() { return minyTot; }
		double getMaxYTotal() { return maxyTot; }
	private:
		double minx,miny,maxx,maxy,sizez,shear;

		// These values are for the whole simulation, not just on this process.
		double minxTot,maxxTot,minyTot,maxyTot;
		ProcessorCommunication &pc;
		SpacingMeasure &spaceMeasure;

		// These are the vector of particle data to send to neighbors.
		// They are encoded as numMoved,numBuffer,x1,y1,z1,vx1,vy1,vz1,x2,...
		std::vector<double> commBuffer[2];
		std::vector<int> removed;
		std::vector<double> sizeBuffer;
};

#endif

#endif
