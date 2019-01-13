/**
 * This header file includes the code to do a forcing based on a tree in parallel.
**/

#include <vector>
#include <algorithm>
#include "ProcessorCommunication.h"

using std::vector;
using std::min;
using std::max;

const int MAX_PARTS=5;

struct KDNode {
	int splitDim;
	double splitVal;
	int firstChild;
	int parts[MAX_PARTS];
	int numParts;
	
	double cmx,cmy,cmz,mass;
	double minx,maxx,miny,maxy,size;
	
	int childNode(double x,double y) {
		if((splitDim==0 && x<splitVal) || (splitDim==1 && y<splitVal)) return firstChild;
		return firstChild+1;
	}
	
	bool operator<(const KDNode &n) { return miny<n.miny; }
};

struct RemoteParticle {
	RemoteParticle() {}
	RemoteParticle(double nx,double ny,double nz,double nmass):x(nx),y(ny),z(nz),mass(nmass) {}
	double x,y,z,mass;
};

class NodeComp {
	public:
		NodeComp(vector<KDNode> &p):pool(p) {}
		bool operator() (int a,int b) {
			return pool[a]<pool[b];
		}
	private:
		vector<KDNode> &pool;
};

template<class Boundary>
class KDGravTree {
	public:
		KDGravTree(double t,ProcessorCommunication &procComm,Boundary &b,double md=0):theta(t),minDist(md),pc(procComm),bounds(b),boundsBuffer(4),sizeBuffer(1) {}
	
		template<class Population>
		void applyForce(Population &pop) {
			pool.resize(2*pop.getNumBodies());
			pool[0].numParts=0;
			firstFree=1;
			for(int i=0; i<pop.getNumReal(); ++i) {
				addParticle(0,pop,i);
			}
			finalize(0,pop);
			printf("Tree finalized. %d\n",pc.getProcessNum());
			communicateTrees(pop);
			printf("Communication done. %d\n",pc.getProcessNum());
			int nb=pop.getNumBodies();
			#pragma omp parallel for schedule(dynamic,100)
			for(int i=0; i<nb; ++i) {
				for(unsigned int j=0; j<receivedTrees.size(); ++j) {
//					printf("%d Doing force on tree %d - %d\n",pc.getProcessNum(),j,receivedTrees[j]);
					doForce(receivedTrees[j],pop,i);
				}
				pop.adjustAfterForce(i);
			}
			printf("Leaving forcing. %d\n",pc.getProcessNum());
		}
	private:
		template<class Population>
		void addParticle(int n,Population &pop,int i) {
			if(pool[n].numParts<0) {
				addParticle(pool[n].childNode(pop.getx(i),pop.gety(i)),pop,i);
			} else if(pool[n].numParts<MAX_PARTS) {
				pool[n].parts[pool[n].numParts]=i;
				++pool[n].numParts;
			} else {
				if(firstFree+1>=pool.size()) {
					pool.resize(pool.size()+pool.size()/2);
				}
				pool[n].firstChild=firstFree;
				pool[firstFree].numParts=0;
				pool[firstFree+1].numParts=0;
				double minx=pop.getx(pool[n].parts[0]);
				double maxx=minx;
				double miny=pop.gety(pool[n].parts[0]);
				double maxy=miny;
				for(int j=1; j<MAX_PARTS; ++j) {
					double x=pop.getx(pool[n].parts[j]);
					if(x<minx) minx=x;
					if(x>maxx) maxx=x;
					double y=pop.gety(pool[n].parts[j]);
					if(y<miny) miny=y;
					if(y>maxy) maxy=y;					
				}
				if(maxx-minx>maxy-miny) {
					pool[n].splitDim=0;
					pool[n].splitVal=0.5*(maxx+minx);
				} else {
					pool[n].splitDim=1;
					pool[n].splitVal=0.5*(maxy+miny);
				}
//				printf("(%e %e) (%e %e) %d %e\n",minx,maxx,miny,maxy,pool[n].splitDim,pool[n].splitVal);
				firstFree+=2;
				pool[n].numParts=-1;
				for(int j=0; j<MAX_PARTS; ++j) {
					addParticle(n,pop,pool[n].parts[j]);
				}
				addParticle(n,pop,i);
			}
		}
	
		template<class Population>
		void finalize(int n,Population &pop) {
			pool[n].mass=0.0;
			pool[n].cmx=0.0;
			pool[n].cmy=0.0;
			pool[n].cmz=0.0;
			if(pool[n].numParts<0) {
				finalize(pool[n].firstChild,pop);
				finalize(pool[n].firstChild+1,pop);
				
				pool[n].mass+=pool[pool[n].firstChild].mass;
				pool[n].cmx+=pool[pool[n].firstChild].mass*pool[pool[n].firstChild].cmx;
				pool[n].cmy+=pool[pool[n].firstChild].mass*pool[pool[n].firstChild].cmy;
				pool[n].cmz+=pool[pool[n].firstChild].mass*pool[pool[n].firstChild].cmz;
				
				pool[n].mass+=pool[pool[n].firstChild+1].mass;
				pool[n].cmx+=pool[pool[n].firstChild+1].mass*pool[pool[n].firstChild+1].cmx;
				pool[n].cmy+=pool[pool[n].firstChild+1].mass*pool[pool[n].firstChild+1].cmy;
				pool[n].cmz+=pool[pool[n].firstChild+1].mass*pool[pool[n].firstChild+1].cmz;

				pool[n].minx=min(pool[pool[n].firstChild].minx,pool[pool[n].firstChild+1].minx);
				pool[n].maxx=max(pool[pool[n].firstChild].maxx,pool[pool[n].firstChild+1].maxx);
				pool[n].miny=min(pool[pool[n].firstChild].miny,pool[pool[n].firstChild+1].miny);
				pool[n].maxy=max(pool[pool[n].firstChild].maxy,pool[pool[n].firstChild+1].maxy);
			} else if(pool[n].numParts>0) {
				pool[n].minx=pop.getx(pool[n].parts[0]);
				pool[n].maxx=pop.getx(pool[n].parts[0]);
				pool[n].miny=pop.gety(pool[n].parts[0]);
				pool[n].maxy=pop.gety(pool[n].parts[0]);
				for(int j=0; j<pool[n].numParts; ++j) {
					double x=pop.getx(pool[n].parts[j]);
					double y=pop.gety(pool[n].parts[j]);
					double mass=pop.getMass(pool[n].parts[j]);
					pool[n].mass+=mass;
					pool[n].cmx+=mass*x;
					pool[n].cmy+=mass*y;
					pool[n].cmz+=mass*pop.getz(pool[n].parts[j]);
					if(x<pool[n].minx) pool[n].minx=x;
					if(x>pool[n].maxx) pool[n].maxx=x;
					if(y<pool[n].miny) pool[n].miny=y;
					if(y>pool[n].maxy) pool[n].maxy=y;
				}
			}
			if(pool[n].mass>0.0) {
				pool[n].cmx/=pool[n].mass;
				pool[n].cmy/=pool[n].mass;
				pool[n].cmz/=pool[n].mass;
				double dx=pool[n].maxx-pool[n].minx;
				double dy=pool[n].maxy-pool[n].miny;
				pool[n].size=max(dx,dy);
			}
		}

		template<class Population>
		void doForce(int n,Population &pop,int i) {
			if(n<0) {
				printf("Negative node in force. %d %d\n",pc.getProcessNum(),n);
				exit(-1);
			}
			if(pool[n].mass<=0.0) return;
//			printf("%d %d %d\n",pc.getProcessNum(),n,i);
			if(pool[n].numParts>0) {
				for(int j=0; j<pool[n].numParts; ++j) {
					int oi=pool[n].parts[j];
					if(oi!=i && oi>=0) {
						double dx=pop.getx(oi)-pop.getx(i);
						double dy=pop.gety(oi)-pop.gety(i);
						double dz=pop.getz(oi)-pop.getz(i);
						double dist=sqrt(dx*dx+dy*dy+dz*dz);
						double mag=pop.getTimeStep()*pop.getMass(oi)/(dist*dist*dist);
						pop.setvx(i,pop.getvx(i)+dx*mag);
						pop.setvy(i,pop.getvy(i)+dy*mag);
						pop.setvz(i,pop.getvz(i)+dz*mag);
					} else if(oi<0) {
//						printf("rp force %d of %d, %e %e %e %e\n",j,pool[n].numParts,remoteParticles[-oi].x,remoteParticles[-oi].y,remoteParticles[-oi].z,remoteParticles[-oi].mass);
						double dx=remoteParticles[-oi].x-pop.getx(i);
						double dy=remoteParticles[-oi].y-pop.gety(i);
						double dz=remoteParticles[-oi].z-pop.getz(i);
						double dist=sqrt(dx*dx+dy*dy+dz*dz);
						if(dist>2*pop.getRadius(i)) {
							double mag=pop.getTimeStep()*remoteParticles[-oi].mass/(dist*dist*dist);
							pop.setvx(i,pop.getvx(i)+dx*mag);
							pop.setvy(i,pop.getvy(i)+dy*mag);
							pop.setvz(i,pop.getvz(i)+dz*mag);
						}
					}
				}
			} else {
				double dx=pool[n].cmx-pop.getx(i);
				double dy=pool[n].cmy-pop.gety(i);
				double dz=pool[n].cmz-pop.getz(i);
				double dsqr=dx*dx+dy*dy+dz*dz;
				if(theta*theta*dsqr>pool[n].size*pool[n].size || pool[n].firstChild<0) {
					double dist=sqrt(dsqr);
					double mag=pop.getTimeStep()*pool[n].mass/(dist*dist*dist);
					pop.setvx(i,pop.getvx(i)+dx*mag);
					pop.setvy(i,pop.getvy(i)+dy*mag);
					pop.setvz(i,pop.getvz(i)+dz*mag);
				} else {
					if(pool[n].firstChild<0) {
						printf("Node %d has children %d\n",n,pool[n].firstChild);
						printf("Distances %e %e\n",theta*theta*dsqr,pool[n].size*pool[n].size);
						return;
					}
					if(dsqr+pool[n].size*pool[n].size<minDist*minDist) return;
					doForce(pool[n].firstChild,pop,i);
					doForce(pool[n].firstChild+1,pop,i);
				}
			}
		}
		
		template<class Population>
		void communicateTrees(Population &pop) {
			receivedTrees.resize(0);
			receivedTrees.push_back(0);
			remoteParticles.resize(0);
			vector<int> rt;
			for(int i=0; i<pc.getNumProcesses(); ++i) {
				int numRead;
				if(i==pc.getProcessNum()) {
					// receive from others.
					boundsBuffer[0]=pool[0].minx;
					boundsBuffer[1]=pool[0].maxx;
					boundsBuffer[2]=pool[0].miny;
					boundsBuffer[3]=pool[0].maxy;
					for(int j=0; j<pc.getNumProcesses(); ++j) {
						if(j!=i) {
							printf("%d recieve from %d\n",i,j);
							pc.sendTo(j,boundsBuffer);
							pc.readFrom(j,sizeBuffer,numRead);
							treeBuffer.resize((int)sizeBuffer[0]);
							pc.readFrom(j,treeBuffer,numRead);
							rt.resize(0);
//							printf("%d reconstruct\n",i);
							reconstructTree(rt);
							for(unsigned int k=0; k<rt.size(); ++k) {
								receivedTrees.push_back(rt[k]);
							}
//							printf("%d going to next\n",i);
						}
					}
				} else {
					// send to others.
					printf("%d send to %d\n",pc.getProcessNum(),i);
					pc.readFrom(i,boundsBuffer,numRead);
//					printf("%d pack tree\n",pc.getProcessNum());
					treeBuffer.resize(0);
					packTree(pop,boundsBuffer[0],boundsBuffer[1],boundsBuffer[2],boundsBuffer[3]);
//					printf("%d done packing tree\n",pc.getProcessNum());
#ifdef AZIMUTHAL_MIRRORS
					double cellWidth=bounds.getMaxYTotal()-bounds.getMinYTotal();
					packTree(pop,boundsBuffer[0],boundsBuffer[1],boundsBuffer[2],boundsBuffer[3],0,cellWidth);
					packTree(pop,boundsBuffer[0],boundsBuffer[1],boundsBuffer[2],boundsBuffer[3],0,-cellWidth);
#endif
					sizeBuffer[0]=treeBuffer.size();
					pc.sendTo(i,sizeBuffer);
					pc.sendTo(i,treeBuffer);
				}
			}
		}
		
		/*
		 * These functions pack and unpack the partial trees.  The format of
		 * the packing is as follows.
		 * 0 - type (0 for CM node, 1 for particle node)
		 * For a CM node you have the following
		 * 1 - cmx
		 * 2 - cmy
		 * 3 - cmz
		 * 4 - cmass
		 * 5 - minx
		 * 6 - maxx
		 * 7 - miny
		 * 8 - maxy
		 * 9 - splitDim
		 * 10 - spitVal
		 * 
		 * For a particle you have the following
		 * 1- numParts
		 * 2 - x
		 * 3 - y
		 * 4 - z
		 * 5 - mass
		 * 6 - x
		 * ...
		 */
		
		template<class Population>
		void packTree(Population &pop,double xmin,double xmax,double ymin,double ymax,double offsetX=0.0,double offsetY=0.0) {
			packNode(pop,offsetX,offsetY,0);
			packRecur(pop,xmin,xmax,ymin,ymax,offsetX,offsetY,0);
		}
		
		template<class Population>
		void packRecur(Population &pop,double xmin,double xmax,double ymin,double ymax,double offsetX,double offsetY,int n) {
//			printf("%d pack recur on %d\n",pc.getProcessNum(),n);
			if(pool[n].numParts<1) {
				// Now check distance and see if we recurse.
				double xSep=std::max(xmin-(pool[n].maxx+offsetX),(pool[n].minx+offsetX)-xmax);
				if(xSep<0.0) xSep=0.0;
				double ySep=std::max(ymin-(pool[n].maxy+offsetY),(pool[n].miny+offsetY)-ymax);
				if(ySep<0.0) ySep=0.0;
				if((xSep*xSep+ySep*ySep)*theta*theta<=pool[n].size*pool[n].size) {
					packNode(pop,offsetX,offsetY,pool[n].firstChild);
					packNode(pop,offsetX,offsetY,pool[n].firstChild+1);
					packRecur(pop,xmin,xmax,ymin,ymax,offsetX,offsetY,pool[n].firstChild);
					packRecur(pop,xmin,xmax,ymin,ymax,offsetX,offsetY,pool[n].firstChild+1);
				}
			}
		}		
		
		template<class Population>
		void packNode(Population &pop,double offsetX,double offsetY,int n) {
			if(pool[n].numParts<1) {
				treeBuffer.push_back(0);
				treeBuffer.push_back(pool[n].cmx+offsetX);
				treeBuffer.push_back(pool[n].cmy+offsetY);
				treeBuffer.push_back(pool[n].cmz);
				treeBuffer.push_back(pool[n].mass);
				treeBuffer.push_back(pool[n].minx+offsetX);
				treeBuffer.push_back(pool[n].maxx+offsetX);
				treeBuffer.push_back(pool[n].miny+offsetY);
				treeBuffer.push_back(pool[n].maxy+offsetY);
				treeBuffer.push_back(pool[n].splitDim);
				double splitOff[]={offsetX,offsetY,0.0};
				treeBuffer.push_back(pool[n].splitVal+splitOff[pool[n].splitDim]);
			} else {
				if(pool[n].numParts==0) {
					printf("There are no particles in %d\n",n);
				}
				treeBuffer.push_back(1);
				treeBuffer.push_back(pool[n].numParts);
				for(int i=0; i<pool[n].numParts; ++i) {
					int p=pool[n].parts[i];
					treeBuffer.push_back(pop.getx(p)+offsetX);
					treeBuffer.push_back(pop.gety(p)+offsetY);
					treeBuffer.push_back(pop.getz(p));
					treeBuffer.push_back(pop.getMass(p));
				}
			}
		}		
		
		void reconstructTree(vector<int> &roots) {
			int curRoot=firstFree;
			roots.push_back(curRoot);
			int addCnt=1;
			for(unsigned int i=0; i<treeBuffer.size();) {
//				printf("%d reconstruct %d of %d\n",pc.getProcessNum(),i,treeBuffer.size());
				if(firstFree>=pool.size()) {
					pool.resize(pool.size()+pool.size()/2);
				}
				pool[firstFree].firstChild=-1;
				if(treeBuffer[i]==0) {
//					printf("CM node\n");
					pool[firstFree].cmx=treeBuffer[i+1];
					pool[firstFree].cmy=treeBuffer[i+2];
					pool[firstFree].cmz=treeBuffer[i+3];
					pool[firstFree].mass=treeBuffer[i+4];
					pool[firstFree].minx=treeBuffer[i+5];
					pool[firstFree].maxx=treeBuffer[i+6];
					pool[firstFree].miny=treeBuffer[i+7];
					pool[firstFree].maxy=treeBuffer[i+8];
					pool[firstFree].splitDim=(int)treeBuffer[i+9];
					pool[firstFree].splitVal=treeBuffer[i+10];
					pool[firstFree].numParts=-1;
					i+=11;
				} else if(treeBuffer[i]==1) {
					pool[firstFree].numParts=(int)treeBuffer[i+1];
//					printf("part node %d\n",pool[firstFree].numParts);
					i+=2;
					for(int j=0; j<pool[firstFree].numParts; ++j) {
						double x=treeBuffer[i];
						double y=treeBuffer[i+1];
						double z=treeBuffer[i+2];
						double mass=treeBuffer[i+3];
						i+=4;
						pool[firstFree].parts[j]=-remoteParticles.size();
						remoteParticles.push_back(RemoteParticle(x,y,z,mass));
//						printf("rp %d of %d, %e %e %e %e\n",j,pool[firstFree].numParts,x,y,z,mass);
						if(j==0) {
							pool[firstFree].minx=x;
							pool[firstFree].maxx=x;
							pool[firstFree].miny=y;
							pool[firstFree].maxy=y;
						} else {
							if(x<pool[firstFree].minx) pool[firstFree].minx=x;
							if(x>pool[firstFree].maxx) pool[firstFree].maxx=x;
							if(y<pool[firstFree].miny) pool[firstFree].miny=y;
							if(y>pool[firstFree].maxy) pool[firstFree].maxy=y;
						}
					}
				} else {
					printf("Error: Bad code in treeBuffer! %d is %f\n",i,treeBuffer[i]);
					exit(-1);
				}
//				printf("Bounds %e %e %e %e\n",
//					pool[firstFree].minx,pool[firstFree].maxx,
//					pool[firstFree].miny,pool[firstFree].maxy);
				pool[firstFree].size=std::max(pool[firstFree].maxx-pool[firstFree].minx,pool[firstFree].maxy-pool[firstFree].miny);
				if(firstFree!=curRoot) {
					if(fitsIn(firstFree,curRoot)) {
						if(addCnt==1) {
//							printf("Add node %d %d\n",curRoot,firstFree);
							addNode(curRoot,firstFree);
							addCnt=0;
//							printf("Done adding\n");
						} else {
							addCnt=1;
						}
					} else {
						if(addCnt==0) {
							printf("Error: new root when expecting second child.\n");
							exit(-1);
						}
//						printf("%d Diff tree %d %d (%e %e) (%e %e) (%e %e) (%e %e)\n",pc.getProcessNum(),curRoot,firstFree,
							pool[firstFree].minx,pool[curRoot].minx,
							pool[firstFree].maxx,pool[curRoot].maxx,
							pool[firstFree].miny,pool[curRoot].miny,
							pool[firstFree].maxy,pool[curRoot].maxy;
						curRoot=firstFree;
						roots.push_back(curRoot);
						addCnt=1;
					}
				} else {
//					printf("%d Diff tree %d %e %e %e %e\n",pc.getProcessNum(),firstFree,
						pool[firstFree].minx,
						pool[firstFree].maxx,
						pool[firstFree].miny,
						pool[firstFree].maxy;
				}
				firstFree++;
			}
		}
		
		/**
		 * Checks if node n fits inside the bounds of node m.
		 */
		bool fitsIn(int n,int m) {
			return pool[n].minx>=pool[m].minx &&
				pool[n].maxx<=pool[m].maxx &&
				pool[n].miny>=pool[m].miny &&
				pool[n].maxy<=pool[m].maxy;
		}
		
		void addNode(int n,int newNode) {
			int p=-1;
			double cx=(pool[newNode].maxx+pool[newNode].minx)*0.5;
			double cy=(pool[newNode].maxy+pool[newNode].miny)*0.5;
			while(n!=-1) {
//				printf("%d %d %d\n",pc.getProcessNum(),n,p);
				p=n;
				n=pool[n].childNode(cx,cy);
				if(n==0) {
					printf("Error %d: add node got to 0, p=%d, newNode=%d, pool.size=%d.\n",pc.getProcessNum(),p,newNode,pool.size());
					printf("p.numParts=%d, n.numParts=%d\n",pool[p].numParts,pool[newNode].numParts);
					printf("cx=%e, cy=%e\n",cx,cy);
					printf("p.splitDim=%d, p.splitVal=%e\n",pool[p].splitDim,pool[p].splitVal);
					printf("p bounds %e %e %e %e\n",pool[p].minx,pool[p].maxx,pool[p].miny,pool[p].maxy);
					printf("new bounds %e %e %e %e\n",pool[newNode].minx,pool[newNode].maxx,pool[newNode].miny,pool[newNode].maxy);
					fflush(stdout);
					exit(-1);
				}
			}
			pool[p].firstChild=newNode;
			if(!fitsIn(newNode,p)) {
				printf("Put %d in %d but it doesn't fit.\n",newNode,p);
					printf("p bounds %e %e %e %e\n",pool[p].minx,pool[p].maxx,pool[p].miny,pool[p].maxy);
					printf("new bounds %e %e %e %e\n",pool[newNode].minx,pool[newNode].maxx,pool[newNode].miny,pool[newNode].maxy);
			}
		}
		
		vector<KDNode> pool;
		int firstFree;
		int finalRoot;
		double theta;
		double minDist;
		int parallelRoot;
		ProcessorCommunication &pc;
		Boundary &bounds;
		vector<double> boundsBuffer;
		vector<double> sizeBuffer;
		vector<double> treeBuffer;
		vector<int> receivedTrees;
		vector<RemoteParticle> remoteParticles;
};

