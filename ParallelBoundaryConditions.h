// ParallelBoundaryConditions.h
// This file contains some classes for different boundary conditions.  These
// classes are one of the main places where code for multiple processors shows
// up, because without gravity using boundaing layers the collisions run like
// they are on a single machine for each timestep.  At the end of the timestep
// the talk to one another about which particles have moved across the lines.
// A significant question here is whether I should make separate versions for
// the sinlge processor BCs or if that should be special code in a more
// complete version.

#ifndef PARALLEL_BOUNDARY_CONDITIONS
#define PARALLEL_BOUNDARY_CONDITIONS

#ifndef BUFFER_WIDTH_MULTIPLE
#define BUFFER_WIDTH_MULTIPLE 3
#endif

#include <vector>
#include "GCPopulation.h"
#include "ProcessorCommunication.h"
#include "Coordinates.h"

#ifdef PARALLEL

class ZeroMinSpacing {
	public:
		double getMinSpacing() {
			return 0.0;
		}
};

template<class SpecificBounds,class SpacingMeasure>
class ParallelBoundary {
	public:
		ParallelBoundary(SpecificBounds &sb,double minX,double maxX,double minY,double maxY,ProcessorCommunication &procComm,SpacingMeasure &sm):specificBounds(sb),minxTot(minX),maxxTot(maxX),minyTot(minY),maxyTot(maxY),pc(procComm),spaceMeasure(sm) {
			minx=minxTot;
			maxx=maxxTot;
			double size=(maxyTot-minyTot)/pc.getNumProcesses();
			miny=minyTot+pc.getProcessNum()*size;
			maxy=minyTot+(pc.getProcessNum()+1)*size;
			timeForLastStep = -1;
			stepCount = 0;
		}

		template<class Population>
		void apply(Population &pop) {
			double bufferSize=BUFFER_WIDTH_MULTIPLE*spaceMeasure.getMinSpacing();
			specificBounds.alterBounds(pop,minx,maxx,miny,maxy,minyTot,maxyTot);
			int numMoved[2]={0,0};
			int numBuffer[2]={0,0};
			commBuffer[0].resize(2);
			commBuffer[1].resize(2);
			removed.resize(0);
			sizeBuffer.resize(2);

			if (timeForLastStep > -1) {
				timeForLastStep = MPI_Wtime() - timeForLastStep;
				for (int i = 0; i < pc.getNumProcesses(); i++) {
//					if (i == pc.getProcessNum()) {
//						fprintf(stderr, "\t%f", timeForLastStep);
//						fflush(stderr);
//					}
					MPI_Barrier(MPI_COMM_WORLD);
				}
				MPI_Barrier(MPI_COMM_WORLD);
//				if (pc.getProcessNum() == 0) {
//					fprintf(stderr, "\n");
//					fflush(stderr);
//				}
				MPI_Barrier(MPI_COMM_WORLD);
			}
#ifdef LOADBALANCE
			if (timeForLastStep > -1) {
				int numRead;
				std::vector<double> loadInfo(2);
				loadInfo[0] = timeForLastStep;
				loadInfo[1] = (miny + maxy) / 2.0;

				//Load-balancing method comes from paper by Leezer and Lewis
				double mid, neighborMid, ratio, newRatio;
				mid = loadInfo[1];
				if (pc.getProcessNum() > 0)
					pc.sendToNeighbor(3, loadInfo);
				if (pc.getProcessNum() < pc.getNumProcesses() - 1) {
					pc.sendToNeighbor(5, loadInfo);
					pc.readFromNeighbor(5, loadInfo, numRead);
					neighborMid = loadInfo[1];
					ratio = (maxy - mid) / (neighborMid - mid);
					newRatio = loadInfo[0]*ratio / (timeForLastStep*(1-ratio) + loadInfo[0]*ratio);
					maxy = newRatio*(neighborMid - mid) + mid;
				}
				if (pc.getProcessNum() > 0) {
					pc.readFromNeighbor(3, loadInfo, numRead);
					neighborMid = loadInfo[1];
					ratio = (miny - neighborMid) / (mid - neighborMid);
					newRatio = timeForLastStep*ratio / (loadInfo[0]*(1-ratio) + timeForLastStep*ratio);
					miny = newRatio*(mid - neighborMid) + neighborMid;
				}
			}
#endif
			bool wrapAround=specificBounds.doWrapAround(minx,maxx,miny,maxy,minyTot,maxyTot);

//			printf("%d - Applying BC %e %e %e %d\n",pc.getProcessNum(),miny,maxy,bufferSize,pop.getNumReal());

			specificBounds.doNonPassedBounds(pop,minx,maxx,miny,maxy,minyTot,maxyTot);

			// Only go up to numReal.  Others are thrown away.
			for(int i=0; i<pop.getNumReal(); ++i) {
				// Check downstream.
//				if(i%10000==0) printf("i[%d]=%e\n",i,pop.gety(i));
				if(pop.gety(i)<miny) {
					if(pc.getProcessNum()==0 && wrapAround) {
						specificBounds.belowMiny(pop,i,minx,maxx,miny,maxy,minyTot,maxyTot);
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
						specificBounds.aboveMaxy(pop,i,minx,maxx,miny,maxy,minyTot,maxyTot);
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
//				printf("%d - Send right\n",pc.getProcessNum());
				sizeBuffer[0]=numMoved[0];
				sizeBuffer[1]=numBuffer[0];
				pc.sendToNeighbor(3,sizeBuffer);
				pc.sendToNeighbor(3,commBuffer[0]);
//				printf("%d - Read left\n",pc.getProcessNum());
				error[0]=pc.readFromNeighbor(5,sizeBuffer,numRead[0]);
				unsigned int size=2+7*(int)(sizeBuffer[0]+sizeBuffer[1]);
				if(size>commBuffer[0].size())
					commBuffer[0].resize(size);
				error[0]=pc.readFromNeighbor(5,commBuffer[0],numRead[0]);
			} else {
//				printf("%d - Read left\n",pc.getProcessNum());
//				printf("%d - Read size\n",pc.getProcessNum());
				error[0]=pc.readFromNeighbor(5,sizeBuffer,numRead[0]);
				unsigned int size=2+7*(int)(sizeBuffer[0]+sizeBuffer[1]);
//				printf("Resize to %d\n",size);
				if(size>recBuffer.size())
					recBuffer.resize(size);
//				printf("%d - Read bodies\n",pc.getProcessNum());
				error[0]=pc.readFromNeighbor(5,recBuffer,numRead[0]);
//				printf("%d - Send right\n",pc.getProcessNum());
				sizeBuffer[0]=numMoved[0];
				sizeBuffer[1]=numBuffer[0];
				pc.sendToNeighbor(3,sizeBuffer);
				pc.sendToNeighbor(3,commBuffer[0]);
				commBuffer[0]=recBuffer;
			}

			// Send left, read right.
			if(pc.getProcessNum()==0) {
//				printf("%d - Send left\n",pc.getProcessNum());
				sizeBuffer[0]=numMoved[1];
				sizeBuffer[1]=numBuffer[1];
				pc.sendToNeighbor(5,sizeBuffer);
				pc.sendToNeighbor(5,commBuffer[1]);
//				printf("%d - Read right\n",pc.getProcessNum());
				error[1]=pc.readFromNeighbor(3,sizeBuffer,numRead[1]);
				unsigned int size=2+7*(int)(sizeBuffer[0]+sizeBuffer[1]);
				if(size>commBuffer[1].size())
					commBuffer[1].resize(size);
				error[1]=pc.readFromNeighbor(3,commBuffer[1],numRead[1]);
			} else {
//				printf("%d - Read right\n",pc.getProcessNum());
				error[1]=pc.readFromNeighbor(3,sizeBuffer,numRead[1]);
				unsigned int size=2+7*(int)(sizeBuffer[0]+sizeBuffer[1]);
				if(size>recBuffer.size())
					recBuffer.resize(size);
				error[1]=pc.readFromNeighbor(3,recBuffer,numRead[1]);
//				printf("%d - Send left\n",pc.getProcessNum());
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
						numBodies++;
					}
				}
			}

			// Cover extra holes if any.
			if(copiedOver<removed.size()) {
				for(unsigned int i=copiedOver; i<removed.size() && removed[i]<=numBodies-1; ++i) {
					bool flag=true;
					while(flag && numBodies-1>=removed[i]) {
						flag=false;
						for(unsigned int j=i; j<removed.size() && !flag; ++j)
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
					numBodies++;
				}
			}
			timeForLastStep = MPI_Wtime();
			stepCount++;
//			printf("%d - Done with BC  %d, %d\n",pc.getProcessNum(),pop.getNumReal(),pop.getNumBodies());
		}

		bool checkIfIn(double x,double y) {
			return specificBounds.checkIfIn(x,y,minx,maxx,miny,maxy,minyTot,maxyTot);
		}

		double getArea() {
			return specificBounds.getArea(minx,maxx,miny,maxy,minyTot,maxyTot);
		}

		double getMinX() { return minx; }
		double getMaxX() { return maxx; }
		double getMinY() { return miny; }
		double getMaxY() { return maxy; }
		double getMinZ() { return specificBounds.getMinZ(); }
		double getMaxZ() { return specificBounds.getMaxZ(); }
		double getMinYTotal() { return minyTot; }
		double getMaxYTotal() { return maxyTot; }
		double getShearOffset() { return specificBounds.getShearOffset();}
	private:
		double minx,miny,maxx,maxy;
		double timeForLastStep;
		int stepCount;

		// These values are for the whole simulation, not just on this process.
		SpecificBounds &specificBounds;
		double minxTot,maxxTot,minyTot,maxyTot;
		ProcessorCommunication &pc;
		SpacingMeasure &spaceMeasure;

		// These are the vector of particle data to send to neighbors.
		// They are encoded as numMoved,numBuffer,x1,y1,z1,vx1,vy1,vz1,x2,...
		std::vector<double> commBuffer[2];
		std::vector<double> recBuffer;
		std::vector<int> removed;
		std::vector<double> sizeBuffer;
};


class SpecificPeriodicWithPhiShift {
	public:
		bool checkIfIn(double x,double y,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return (x>=minx && x<=maxx && y>=miny && y<=maxy);
		}

		double getArea(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return (maxx-minx)*(maxy-miny);
		}

		template<class Population>
		void alterBounds(Population &pop,double &minx,double &maxx,double &miny,double &maxy,double &minyTot,double &maxyTot) {
			miny-=2.0*GCCoords::A0*minx*pop.getTimeStep();
			maxy-=2.0*GCCoords::A0*minx*pop.getTimeStep();
			minyTot-=2.0*GCCoords::A0*minx*pop.getTimeStep();
			maxyTot-=2.0*GCCoords::A0*minx*pop.getTimeStep();
		}

		bool doWrapAround(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			double deltay=4.0*3.141592654*GCCoords::A0*maxx;
			double tmpy=miny;
			while(tmpy<-3.141592654) tmpy+=6.283185307;
			return tmpy<-2.0*deltay || tmpy>deltay;
		}

		template<class Population>
		void doNonPassedBounds(Population &pop,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		template<class Population>
		void belowMiny(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			pop.setY(i,pop.getY(i)+(maxyTot-minyTot));
			pop.setPhi(i,pop.getPhi(i)-(maxyTot-minyTot)/(2.0*GCCoords::A0*pop.getX(i)));
			pop.setCartAfterForce(i);
		}

		template<class Population>
		void aboveMaxy(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			pop.setY(i,pop.getY(i)-(maxyTot-minyTot));
			pop.setPhi(i,pop.getPhi(i)+(maxyTot-minyTot)/(2.0*GCCoords::A0*pop.getX(i)));
			pop.setCartAfterForce(i);
		}

		template<class Population>
		void belowMinx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		template<class Population>
		void aboveMaxx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		double getShearOffset() { return 0.0;}
};


class SpecificSingleOrbitAzimuthal {
	public:
		bool checkIfIn(double x,double y,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return true;
		}

		double getArea(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return 1.0;
		}

		template<class Population>
		void alterBounds(Population &pop,double &minx,double &maxx,double &miny,double &maxy,double &minyTot,double &maxyTot) {
			miny-=2.0*GCCoords::A0*minx*pop.getTimeStep();
			maxy-=2.0*GCCoords::A0*minx*pop.getTimeStep();
			minyTot-=2.0*GCCoords::A0*minx*pop.getTimeStep();
		}

		bool doWrapAround(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			double deltay=4.0*3.141592654*GCCoords::A0*maxx;
			double tmpy=miny;
			while(tmpy<-3.141592654) tmpy+=6.283185307;
			return tmpy<-1.5*deltay || tmpy>deltay*0.5;
		}

		template<class Population>
		void doNonPassedBounds(Population &pop,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		template<class Population>
		void belowMiny(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			double yTotDiff=4.0*3.141592654*GCCoords::A0*pop.getX(i);
			pop.setY(i,pop.getY(i)+yTotDiff);
			pop.setCartAfterForce(i);
		}

		template<class Population>
		void aboveMaxy(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			double yTotDiff=4.0*3.141592654*GCCoords::A0*pop.getX(i);
			pop.setY(i,pop.getY(i)-yTotDiff);
			pop.setCartAfterForce(i);
		}

		template<class Population>
		void belowMinx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		template<class Population>
		void aboveMaxx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		double getShearOffset() { return 0.0;}
};

class SpecificFixedPeriodic {
	public:
		SpecificFixedPeriodic(bool zero=false):zeroeOnWrap(zero) {
		}

		bool checkIfIn(double x,double y,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return (x>=minx && x<=maxx && y>=miny && y<=maxy);
		}

		double getArea(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return (maxx-minx)*(maxy-miny);
		}

		template<class Population>
		void alterBounds(Population &pop,double &minx,double &maxx,double &miny,double &maxy,double &minyTot,double &maxyTot) {
		}

		bool doWrapAround(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return true;
		}

		template<class Population>
		void doNonPassedBounds(Population &pop,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		template<class Population>
		void belowMiny(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			pop.setY(i,pop.getY(i)+(maxyTot-minyTot));
			if(zeroeOnWrap) pop.sete(i,0);
			pop.setCartAfterForce(i);
		}

		template<class Population>
		void aboveMaxy(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			pop.setY(i,pop.getY(i)-(maxyTot-minyTot));
			pop.setCartAfterForce(i);
		}

		template<class Population>
		void belowMinx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		template<class Population>
		void aboveMaxx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		double getShearOffset() { return 0.0;}
	private:
		bool zeroeOnWrap;
};

class SpecificFixedSingleOrbitAzimuthal {
	public:
		bool checkIfIn(double x,double y,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return true;
		}

		double getArea(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return 1.0;
		}

		template<class Population>
		void alterBounds(Population &pop,double &minx,double &maxx,double &miny,double &maxy,double &minyTot,double &maxyTot) {
		}

		bool doWrapAround(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			double deltay=4.0*3.141592654*GCCoords::A0*maxx;
			double tmpy=miny;
			while(tmpy<-3.141592654) tmpy+=6.283185307;
			return tmpy<-1.5*deltay || tmpy>deltay*0.5;
		}

		template<class Population>
		void doNonPassedBounds(Population &pop,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		template<class Population>
		void belowMiny(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			double yTotDiff=4.0*3.141592654*GCCoords::A0*pop.getX(i);
			pop.setY(i,pop.getY(i)+yTotDiff);
			pop.setCartAfterForce(i);
		}

		template<class Population>
		void aboveMaxy(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			double yTotDiff=4.0*3.141592654*GCCoords::A0*pop.getX(i);
			pop.setY(i,pop.getY(i)-yTotDiff);
			pop.setCartAfterForce(i);
		}

		template<class Population>
		void belowMinx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		template<class Population>
		void aboveMaxx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		double getShearOffset() { return 0.0;}
};

class SpecificShearedPeriodic {
	public:
		SpecificShearedPeriodic(double sz,double sh):sizez(sz),shear(sh) {}

		bool checkIfIn(double x,double y,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return (x>=minx && x<=maxx && y>=miny && y<=maxy);
		}

		double getArea(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return (maxx-minx)*(maxy-miny);
		}

		template<class Population>
		void alterBounds(Population &pop,double &minx,double &maxx,double &miny,double &maxy,double &minyTot,double &maxyTot) {
		}

		bool doWrapAround(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return true;
		}

		template<class Population>
		void doNonPassedBounds(Population &pop,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
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
			}
		}

		template<class Population>
		void belowMiny(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			pop.sety(i,pop.gety(i)+(maxyTot-minyTot));
			pop.setvx(i,pop.getvx(i)+shear);
		}

		template<class Population>
		void aboveMaxy(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			pop.sety(i,pop.gety(i)-(maxyTot-minyTot));
			pop.setvx(i,pop.getvx(i)-shear);
		}

		template<class Population>
		void belowMinx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		template<class Population>
		void aboveMaxx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		double getMinZ() { return -sizez; }
		double getMaxZ() { return sizez; }
		double getShearOffset() { return 0.0;}
	private:
		double sizez,shear;
};

template<class MoonClass>
class SpecificFixedPeriodicWithMoon {
	public:
		SpecificFixedPeriodicWithMoon(bool zero=false):zeroeOnWrap(zero),moon(0) {
		}

		bool checkIfIn(double x,double y,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return (x>=minx && x<=maxx && y>=miny && y<=maxy);
		}

		double getArea(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return (maxx-minx)*(maxy-miny);
		}

		void setMoon(MoonClass *m) {
			moon=m;
		}

		template<class Population>
		void alterBounds(Population &pop,double &minx,double &maxx,double &miny,double &maxy,double &minyTot,double &maxyTot) {
		}

		bool doWrapAround(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			if(moon!=0) {
				return moon->doWrapAround(minyTot,maxyTot);
			}
			return true;
		}

		template<class Population>
		void doNonPassedBounds(Population &pop,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		template<class Population>
		void belowMiny(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			pop.setY(i,pop.getY(i)+(maxyTot-minyTot));
			if(zeroeOnWrap) pop.sete(i,0);
			pop.setCartAfterForce(i);
		}

		template<class Population>
		void aboveMaxy(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			pop.setY(i,pop.getY(i)-(maxyTot-minyTot));
			pop.setCartAfterForce(i);
		}

		template<class Population>
		void belowMinx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		template<class Population>
		void aboveMaxx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		double getShearOffset() { return 0.0;}
	private:
		bool zeroeOnWrap;
		MoonClass *moon;
};

template<class MoonClass>
class SpecificFixedSingleOrbitAzimuthWithMoon {
	public:
		SpecificFixedSingleOrbitAzimuthWithMoon(bool zero=false):zeroeOnWrap(zero),moon(0) {
		}

		bool checkIfIn(double x,double y,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return true;
		}

		double getArea(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return 1.0;
		}

		void setMoon(MoonClass *m) {
			moon=m;
		}

		template<class Population>
		void alterBounds(Population &pop,double &minx,double &maxx,double &miny,double &maxy,double &minyTot,double &maxyTot) {
		}

		bool doWrapAround(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			if(moon!=0) {
				return moon->doWrapAround(minyTot,maxyTot);
			}
			return true;
		}

		template<class Population>
		void doNonPassedBounds(Population &pop,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		template<class Population>
		void belowMiny(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			pop.setY(i,pop.getY(i)+(maxyTot-minyTot));
			if(zeroeOnWrap) pop.sete(i,0);
			pop.setCartAfterForce(i);
		}

		template<class Population>
		void aboveMaxy(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			pop.setY(i,pop.getY(i)-(maxyTot-minyTot));
			pop.setCartAfterForce(i);
		}

		template<class Population>
		void belowMinx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		template<class Population>
		void aboveMaxx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		double getShearOffset() { return 0.0;}
	private:
		bool zeroeOnWrap;
		MoonClass *moon;
};

class SpecificSlidingBrick {
	public:
		SpecificSlidingBrick():shearOffset(0.0) {}

		bool checkIfIn(double x,double y,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return (x>=minx && x<=maxx && y>=miny && y<=maxy);
		}

		double getArea(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return (maxx-minx)*(maxy-miny);
		}

		template<class Population>
		void alterBounds(Population &pop,double &minx,double &maxx,double &miny,double &maxy,double &minyTot,double &maxyTot) {
		}

		bool doWrapAround(double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			return true;
		}

		template<class Population>
		void doNonPassedBounds(Population &pop,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			for(int i=0; i<pop.getNumReal(); ++i) {
				while(pop.getx(i)<minx) {
//					printf("Low side %e\n",pop.getx(i));
					pop.setx(i,pop.getx(i)+(maxx-minx));
//					printf("After %e\n",pop.getx(i));
					double newy=pop.gety(i)+2*shearOffset;
					while(newy>maxyTot) newy-=(maxyTot-minyTot);
					while(newy<minyTot) newy+=(maxyTot-minyTot);
					pop.sety(i,newy);
					pop.setvy(i,pop.getvy(i)-2.0*GCCoords::A0*(maxx-minx));
					pop.adjustAfterForce(i);
				}
				while(pop.getx(i)>maxx) {
//					printf("High side %e\n",pop.getx(i));
					pop.setx(i,pop.getx(i)-(maxx-minx));
//					printf("After %e\n",pop.getx(i));
					double newy=pop.gety(i)-2*shearOffset;
					while(newy>maxyTot) newy-=(maxyTot-minyTot);
					while(newy<minyTot) newy+=(maxyTot-minyTot);
					pop.sety(i,newy);
					pop.setvy(i,pop.getvy(i)+2.0*GCCoords::A0*(maxx-minx));
					pop.adjustAfterForce(i);
				}
			}

            shearOffset -= 2*GCCoords::A0*(maxx-minx)*0.5*pop.getTimeStep();
            if(shearOffset < minyTot) { shearOffset += maxyTot-minyTot; }
		}

		template<class Population>
		void belowMiny(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			pop.setY(i,pop.getY(i)+(maxyTot-minyTot));
			pop.setCartAfterForce(i);
		}

		template<class Population>
		void aboveMaxy(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
			pop.setY(i,pop.getY(i)-(maxyTot-minyTot));
			pop.setCartAfterForce(i);
		}

		template<class Population>
		void belowMinx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		template<class Population>
		void aboveMaxx(Population &pop,int i,double minx,double maxx,double miny,double maxy,double minyTot,double maxyTot) {
		}

		double getShearOffset() { return shearOffset;}

	private:
        double shearOffset;
};



#endif

#endif
