/**
 * This header file includes the code to do a forcing based on a tree in
 * parallel. This is my version of a spatial KD tree used for collision
 * detection and gravity.  Note that the build happens in the collision
 * detection so when a double force is built the collision detection has to
 * come first.
**/

#include <vector>
#include <algorithm>
#include "ProcessorCommunication.h"

#ifndef PARALLEL_GRAV_COLL_TREE
#define PARALLEL_GRAV_COLL_TREE

const int MAX_NUM=3;
const int DIM=3;

using std::vector;
using std::min;
using std::max;

struct KDBounds {
	double min[DIM],max[DIM];
};

struct KDNode {
	int splitDim;
	double splitVal;
	int firstChild;
	int parts[MAX_NUM];
	int numParts;
	
	double cmx,cmy,cmz,mass;
	KDBounds bounds;
	double size;
	double maxRad;
	double cnt;
	double velDispSum;
	double searchRadius;
	
	int childNode(double x,double y,double z) {
		if((splitDim==0 && x<=splitVal) || (splitDim==1 && y<=splitVal) || (splitDim==2 && z<=splitVal)) return firstChild;
		return firstChild+1;
	}
	
	bool operator<(const KDNode &n) { return bounds.min[1]<n.bounds.min[1]; }
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

struct Location {
	double c[DIM];
};

struct RecurseInfo {
	int n1,n2;
	KDBounds b1,b2;
	int end;
};

struct GravRecurseInfo {
	GravRecurseInfo(int n1,int n2,double ox,double oy):puller(n1),pulled(n2),offsetX(ox),offsetY(oy) {}
	int puller,pulled;
	double offsetX,offsetY;
};

struct AccelVect {
	AccelVect():ax(0.0),ay(0.0),az(0.0) {}
	double ax,ay,az;
};

template<class Boundary>
class ParallelGravCollTree {
	public:
		ParallelGravCollTree(double t,ProcessorCommunication &procComm,Boundary &b,double md=0):acc(omp_get_max_threads()),theta(t),pc(procComm),bounds(b),minDist(md),boundsBuffer(6),sizeBuffer(1) {
			maxDepth = (int) (log(omp_get_max_threads())/log(2) + 2);
			omp_init_lock(&lock);
		}

		~ParallelGravCollTree() {
			omp_destroy_lock(&lock);
		}
	
		template<class Population>
		void applyForce(Population &pop) {
			for(int i=0; i<acc.size(); ++i) {
				acc[i].resize(pop.getNumBodies());
			}
			printf("Do tree gravity.\n");
			fflush(stdout);
			buildReal(pop);
			double cellWidth=bounds.getMaxYTotal()-bounds.getMinYTotal();
#ifdef FULL_MIRRORS
            double cellHeight=bounds.getMaxX()-bounds.getMinX();
            double shearOffset=bounds.getShearOffset();
#endif
			communicateTrees(pop);
			printf("Communication done. %d\n",pc.getProcessNum());
			fflush(stdout);
			int nb=pop.getNumReal();
#ifdef PFORCE
			int nt=acc.size();
			#pragma omp parallel for
			for(int i=0; i<nt; ++i) {
				for(int j=0; j<nb; ++j) {
					acc[i][j].ax=0.0;
					acc[i][j].ay=0.0;
					acc[i][j].az=0.0;
				}
			}
			vector<GravRecurseInfo> stack;
			doPPForce(0,0,pop,0,0,stack,0,acc[0]);
#ifdef AZIMUTHAL_MIRRORS
			for(int j=1; j<=AZIMUTHAL_MIRRORS; ++j) {
				doPPForce(0,0,pop,0,cellWidth*j,stack,0,acc[0]);
				doPPForce(0,0,pop,0,-cellWidth*j,stack,0,acc[0]);
#ifdef FULL_MIRRORS
				for(int k=1; k<=FULL_MIRRORS; ++k) {
					double totalShear=k*shearOffset;
					while(totalShear<bounds.getMinY()) totalShear+=cellWidth;
					doPPForce(0,0,pop,cellHeight*k,totalShear+cellWidth*j,stack,0,acc[0]);
					doPPForce(0,0,pop,cellHeight*k,totalShear-cellWidth*j,stack,0,acc[0]);
					doPPForce(0,0,pop,-cellHeight*k,-totalShear+cellWidth*j,stack,0,acc[0]);
					doPPForce(0,0,pop,-cellHeight*k,-totalShear-cellWidth*j,stack,0,acc[0]);
				}
#endif
			}
#ifdef FULL_MIRRORS
			for(int k=1; k<=FULL_MIRRORS; ++k) {
				double totalShear=k*shearOffset;
				while(totalShear<bounds.getMinY()) totalShear+=cellWidth;
				doPPForce(0,0,pop,cellHeight*k,totalShear,stack,0,acc[0]);
				doPPForce(0,0,pop,-cellHeight*k,-totalShear,stack,0,acc[0]);
			}
#endif
#endif
			for(unsigned int j=0; j<receivedTrees.size(); ++j) {
				doPPForce(receivedTrees[j],0,pop,0,0,stack,0,acc[0]);
			}
			int ssize = (int) stack.size();
			#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < ssize; ++i) {
				doPPForce(stack[i].puller,stack[i].pulled,pop,stack[i].offsetX,stack[i].offsetY,stack,-1000000,acc[omp_get_thread_num()]);
			}
			#pragma omp parallel for
			for(int j=0; j<nb; ++j) {
				for(int i=0; i<acc.size(); ++i) {
					pop.setvx(j,pop.getvx(j)+acc[i][j].ax);
					pop.setvy(j,pop.getvy(j)+acc[i][j].ay);
					pop.setvz(j,pop.getvz(j)+acc[i][j].az);
				}
				pop.adjustAfterForce(j);
			}
#else
			#pragma omp parallel for schedule(dynamic,100)
			for(int i=0; i<nb; ++i) {
				doForce(0,pop,i,0,0);
#ifdef AZIMUTHAL_MIRRORS
	            doForce(0,pop,i,0,cellWidth);
	            doForce(0,pop,i,0,-cellWidth);
#endif
#ifdef FULL_MIRRORS
	            doForce(0,pop,i,cellHeight,shearOffset+cellWidth);
	            doForce(0,pop,i,cellHeight,shearOffset);
	            doForce(0,pop,i,cellHeight,shearOffset-cellWidth);
	            doForce(0,pop,i,-cellHeight,-shearOffset+cellWidth);
	            doForce(0,pop,i,-cellHeight,-shearOffset);
	            doForce(0,pop,i,-cellHeight,-shearOffset-cellWidth);
#endif
				for(unsigned int j=0; j<receivedTrees.size(); ++j) {
					doForce(receivedTrees[j],pop,i,0,0);
				}
				pop.adjustAfterForce(i);
			}
#endif
		}

		template<class Population>
		void build(Population &pop) {
#ifndef GRAVITY
			buildReal(pop);
#endif
		}
		
		template<class Population>
		void buildReal(Population &pop) {
			printf("Building.\n");
			pool.resize(pop.getNumBodies()*2);
			location.resize(pop.getNumBodies());
			pool[0].parts[0]=0;
			pool[0].numParts=0;
			firstFree=1;
			for(int i=0; i<DIM; ++i) {
				totalBounds.min[i]=1e100;
				totalBounds.max[i]=-1e100;
			}

			for(int i=0; i<pop.getNumBodies(); ++i) {
				for(int j=0; j<DIM; ++j) {
					location[i].c[j] = pop.get(i,j);
					if (pop.get(i,j) < totalBounds.min[j])
						totalBounds.min[j]=pop.get(i,j);
					if (pop.get(i,j) > totalBounds.max[j])
						totalBounds.max[j]=pop.get(i,j);
				}
			}

			vector<int> index(pop.getNumBodies());
			for(int i=0; i<pop.getNumBodies(); ++i) {
				index[i]=i;
			}

			//printf("Build recur to %d.\n",maxDepth);
			//fflush(stdout);
			vector<RecurseInfo> info(0);
			buildRecur(pop,index,0,pop.getNumReal(),0,totalBounds,info);
			int isize = (int) info.size();
			//printf("Build loop. %d\n",isize);
			//fflush(stdout);
			#pragma omp parallel for schedule(dynamic,100)
			for (int i = 0; i < isize; i++)
				buildRecur(pop,index,info[i].n2,info[i].end,info[i].n1,info[i].b1,info,-1000000);

			info.resize(0);
			//printf("finalize top.\n");
			//fflush(stdout);
			finalize(0,pop,info);
			isize = (int) info.size();
			//printf("finalize loop. %d\n",isize);
			//fflush(stdout);
			#pragma omp parallel for schedule(dynamic,100)
			for (int i = 0; i < isize; i++)
				finalize(info[i].n1,pop,info,-1000000);
			info.resize(0);
			//printf("finalize end.\n");
			//fflush(stdout);
			finalize(0,pop,info);

			//validateTree(0,pop);

			printf("End build %d nodes\n",firstFree);
			fflush(stdout);
//			printTree(pop,0);
		}

		double getMinSpacing() {
			return 3*(2*pool[0].maxRad+pool[0].searchRadius);
		}

		template<class Population,class CollisionManager>
		void findAllCollisions(Population &pop,CollisionManager &cm) {
			//printf("Find all.\n");
			KDBounds b1,b2;
			for(int i=0; i<DIM; ++i) {
				b1.min[i]=totalBounds.min[i];
				b1.max[i]=totalBounds.max[i];
				b2.min[i]=totalBounds.min[i];
				b2.max[i]=totalBounds.max[i];
			}		
			buildLockGrid(pop);

			//Recurse through the top maxDepth levels of the tree.
			//Instead of making further recursive calls, store the parameters for
			//each call as a RecurseInfo and store it in the info queue.

			//printf("Top.\n");
			std::vector<RecurseInfo> info(0);
			recurseForAll(pop,cm,0,b1,0,b2,info);
			
			//Make all stored recursive calls in parallel.

			//printf("Others.\n");
			int isize = (int) info.size();
			#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < isize; i++)
				recurseForAll(pop,cm,info[i].n1,info[i].b1,info[i].n2,info[i].b2,info,-1000000);
//			printf("End all\n");
		}

		template<class Population,class CollisionManager>
		void findCollisions(Population &pop,CollisionManager &cm,int p,int ignore) {
			KDBounds b;
			for(int i=0; i<DIM; ++i) {
				b.min[i]=totalBounds.min[i];
				b.max[i]=totalBounds.max[i];
			}			
			recurseForOne(pop,cm,p,ignore,b,0);
		}

		int getBinX(int p) {
			return (int) ((location[p].c[0] - totalBounds.min[0])/gridSpacing);
		}

		int getBinY(int p) {
			return (int) ((location[p].c[1] - totalBounds.min[1])/gridSpacing);
		}

		bool isSafe(int x, int y) {
			int startx = (x>0) ? -1 : 0;
			int endx = (x<(int)inUse.size()-1) ? 2 : 1;
			int starty = (y>0) ? -1 : 0;
			int endy = (y<(int)inUse[0].size()-1) ? 2 : 1;
			for (int i = startx; i < endx; i++)
			{
				for (int j = starty; j < endy; j++)
				{
					if (inUse[x+i][y+j])
						return false;
				}
			}
			return true;
		}

		void setInUse(int x, int y, bool val) {
			inUse[x][y] = val;
		}

	private:
		//To synchronize the KDTree we must create a grid back-end to keep track
		//of which regions contain particles that are currently being processed.
		template<class Population>
		void buildLockGrid(Population &pop) {
			gridSpacing = 2*pool[0].maxRad + pool[0].searchRadius;
			int xlen = (int) ((totalBounds.max[0] - totalBounds.min[0])/gridSpacing + 1);
			int ylen = (int) ((totalBounds.max[1] - totalBounds.min[1])/gridSpacing + 1);
			//printf("%d (%e %e %e) %d (%e %e)\n",xlen,totalBounds.max[0],totalBounds.min[0],gridSpacing,ylen,totalBounds.max[1],totalBounds.min[1]);
			inUse.resize(xlen);
			for (int i = 0; i < xlen; i++)
			{
				inUse[i].resize(ylen);
				for (int j = 0; j < ylen; j++)
					inUse[i][j] = false;
			}
		}

		/*void printLeafInfo(int node) {
			if (pool[node].numParts > 0) {
				double spacing = 2*pool[node].maxRad + pool[node].searchRadius;
				fprintf(stderr, "%d leaf[%d] : gridspacing=%e cnt=%e sum=%e\n", firstMade->getProcessNum(), node, spacing, pool[node].cnt, pool[node].velDispSum);
				fflush(stderr);
			}
			else {
				printLeafInfo(pool[node].firstChild);
				printLeafInfo(pool[node].firstChild + 1);
			}
		}*/

		template<class Population>
		void buildRecur(Population &pop,vector<int> &index,int start,int end,int node,KDBounds &bounds,vector<RecurseInfo> &info,int depth=0) {
			if (depth >= maxDepth)
			{
				info.resize(info.size()+1);
				info[info.size()-1].n1=node;
				info[info.size()-1].n2=start;
				info[info.size()-1].end=end;
				info[info.size()-1].b1=bounds;
				return;
			}
			//printf("br %d %d %d %d\n",start,end,node,depth);
			if(end-start<=MAX_NUM) {
				//printf("%d-%d=%d\n",end,start,end-start);
				pool[node].numParts=end-start;
				for(int i=0; i<pool[node].numParts; ++i) {
					pool[node].parts[i]=index[start+i];
				}
				pool[node].firstChild=-100;
			} else {
				int mid=(start+end)/2;
				int dim=0;
				for(int i=1; i<DIM; ++i) {
					if(bounds.max[i]-bounds.min[i]>bounds.max[dim]-bounds.min[dim]) dim=i;
				}
				//printf("Split %d %d in %d\n",end,start,dim);
				findMid(pop,index,start,end,mid,dim);
			/*	for(int i=start; i<end; ++i) {
					if((i<=mid) ^ (pop.get(index[i],dim)<=pop.get(index[mid],dim))) {
						printf("Error!!! %d %d %d %e %e\n",i,mid,dim,pop.get(index[i],dim),pop.get(index[mid],dim));
					}
				}*/
				pool[node].numParts=-1;
				omp_set_lock(&lock);
				pool[node].firstChild=firstFree;
				firstFree+=2;
				omp_unset_lock(&lock);
				pool[node].splitDim=dim;
				pool[node].splitVal=pop.get(index[mid],dim);
				double tmp=bounds.max[dim];
				bounds.max[dim]=pool[node].splitVal;
				buildRecur(pop,index,start,mid+1,pool[node].firstChild,bounds,info,depth+1);
				bounds.max[dim]=tmp;
				tmp=bounds.min[dim];
				bounds.min[dim]=pool[node].splitVal;
				buildRecur(pop,index,mid+1,end,pool[node].firstChild+1,bounds,info,depth+1);
				bounds.min[dim]=tmp;
			}
		}

		template<class Population>
		void findMid(Population &pop,vector<int> &index,int start,int end,int mid,int dim) {
			while(start<end) {
				//printf("Find %d %d %d\n",start,end,mid);
				int low=start+1;
				int high=end-1;
				double pivot=pop.get(index[start],dim);
				while(low<=high) {
					if(pop.get(index[low],dim)<=pivot) {
						low++;
					} else {
						std::swap(index[low],index[high]);
						high--;
					}
				}
				std::swap(index[start],index[high]);
				if(high==mid) return;
				if(high<mid) {
					start=high+1;
				} else {
					end=high;
				}
			}
		}

		template<class Population>
		void validateTree(int n,Population &pop) {
			if(pool[n].numParts<0) {
				int fc=pool[n].firstChild;
				for(int i=0; i<DIM; ++i) {
					if(pool[fc].bounds.min[i]<pool[n].bounds.min[i])
						printf("Invalid min0 node=%d child=%d dim=%n",n,fc,i);
					if(pool[fc].bounds.max[i]>pool[n].bounds.max[i])
						printf("Invalid min0 node=%d child=%d dim=%n",n,fc,i);
					if(pool[fc+1].bounds.min[i]<pool[n].bounds.min[i])
						printf("Invalid min1 node=%d child=%d dim=%n",n,fc+1,i);
					if(pool[fc+1].bounds.max[i]>pool[n].bounds.max[i])
						printf("Invalid min1 node=%d child=%d dim=%n",n,fc+1,i);
				}
				if(pool[fc].bounds.min[pool[n].splitDim]>pool[n].splitVal)
					printf("Invalid side of split min0 node=%d child=%d dim=%d\n",n,fc,pool[n].splitDim);
				if(pool[fc].bounds.max[pool[n].splitDim]>pool[n].splitVal)
					printf("Invalid side of split max0 node=%d child=%d dim=%d\n",n,fc,pool[n].splitDim);
				if(pool[fc+1].bounds.min[pool[n].splitDim]<pool[n].splitVal)
					printf("Invalid side of split min1 node=%d child=%d dim=%d\n",n,fc+1,pool[n].splitDim);
				if(pool[fc+1].bounds.max[pool[n].splitDim]<pool[n].splitVal)
					printf("Invalid side of split max1 node=%d child=%d dim=%d\n",n,fc+1,pool[n].splitDim);
				if(pool[fc].numParts>0) {
					for(int i=0; i<pool[fc].numParts; ++i) {
						if(pop.get(pool[fc].parts[i],pool[n].splitDim)>pool[n].splitVal)
							printf("Invalid: particle %d above split.\n",pool[fc].parts[i]);
					}
				}
				if(pool[fc+1].numParts>0) {
					for(int i=0; i<pool[fc+1].numParts; ++i) {
						if(pop.get(pool[fc+1].parts[i],pool[n].splitDim)<pool[n].splitVal)
							printf("Invalid: particle %d below split.\n",pool[fc].parts[i]);
					}
				}
				validateTree(fc,pop);
				validateTree(fc+1,pop);
			} else if(pool[n].numParts>0) {
				for(int i=0; i<pool[n].numParts; ++i) {
					if(pop.getx(pool[n].parts[i])<pool[n].bounds.min[0])
						printf("Invalid x node=%d part=%d\n",n,pool[n].parts[i]);
					if(pop.getx(pool[n].parts[i])>pool[n].bounds.max[0])
						printf("Invalid x node=%d part=%d\n",n,pool[n].parts[i]);
					if(pop.gety(pool[n].parts[i])<pool[n].bounds.min[1])
						printf("Invalid y node=%d part=%d\n",n,pool[n].parts[i]);
					if(pop.gety(pool[n].parts[i])>pool[n].bounds.max[1])
						printf("Invalid y node=%d part=%d\n",n,pool[n].parts[i]);
					if(pop.getz(pool[n].parts[i])<pool[n].bounds.min[2])
						printf("Invalid z node=%d part=%d\n",n,pool[n].parts[i]);
					if(pop.getz(pool[n].parts[i])>pool[n].bounds.max[2])
						printf("Invalid z node=%d part=%d\n",n,pool[n].parts[i]);
				}
			}
		}

		template<class Population>
		void finalize(int n,Population &pop,vector<RecurseInfo> &info,int depth=0) {
			if (depth >= maxDepth)
			{
				info.resize(info.size()+1);
				info[info.size()-1].n1=n;
				return;
			}
#ifdef GRAVITY
			pool[n].mass=0.0;
			pool[n].cmx=0.0;
			pool[n].cmy=0.0;
			pool[n].cmz=0.0;
#endif
			if(pool[n].numParts<0) {
				pool[n].velDispSum=0.0;
				pool[n].cnt=0.0;
				for(int j=0; j<2; ++j) {
					finalize(pool[n].firstChild+j,pop,info,depth+1);
					pool[n].velDispSum+=pool[pool[n].firstChild+j].velDispSum;
					pool[n].cnt+=pool[pool[n].firstChild+j].cnt;
#ifdef GRAVITY
					pool[n].mass+=pool[pool[n].firstChild+j].mass;
					pool[n].cmx+=pool[pool[n].firstChild+j].mass*pool[pool[n].firstChild+j].cmx;
					pool[n].cmy+=pool[pool[n].firstChild+j].mass*pool[pool[n].firstChild+j].cmy;
					pool[n].cmz+=pool[pool[n].firstChild+j].mass*pool[pool[n].firstChild+j].cmz;
#endif
				}
				for(int j=0; j<DIM; ++j) {
					pool[n].bounds.min[j]=min(pool[pool[n].firstChild].bounds.min[j],pool[pool[n].firstChild+1].bounds.min[j]);
					pool[n].bounds.max[j]=max(pool[pool[n].firstChild].bounds.max[j],pool[pool[n].firstChild+1].bounds.max[j]);
				}
				pool[n].maxRad=max(pool[pool[n].firstChild].maxRad,pool[pool[n].firstChild+1].maxRad);
			} else if(pool[n].numParts>0) {
				pool[n].maxRad=0.0;
				pool[n].velDispSum=0.0;
				pool[n].cnt=0.0;
				for(int i=0; i<pool[n].numParts; ++i) {
					double radius=pop.getRadius(pool[n].parts[i]);
					if(radius>pool[n].maxRad) pool[n].maxRad=radius;
					for(int j=i+1; j<pool[n].numParts; ++j) {
						double dx=pop.getvx(pool[n].parts[i])-pop.getvx(pool[n].parts[j]);
						double dy=pop.getvy(pool[n].parts[i])-pop.getvy(pool[n].parts[j]);
						double dz=pop.getvz(pool[n].parts[i])-pop.getvz(pool[n].parts[j]);
						pool[n].velDispSum+=dx*dx+dy*dy+dz*dz;
						pool[n].cnt+=1.0;
//						printf("%e %e %e %e %e\n",dx,dy,dz,searchRadius,cnt);
					}
				}
				for(int j=0; j<DIM; ++j) {
					pool[n].bounds.min[j]=pop.get(pool[n].parts[0],j);
					pool[n].bounds.max[j]=pop.get(pool[n].parts[0],j);
				}
				for(int j=0; j<pool[n].numParts; ++j) {
					double x[3]={pop.getx(pool[n].parts[j]),pop.gety(pool[n].parts[j]),pop.getz(pool[n].parts[j])};
#ifdef GRAVITY
					if(pool[n].parts[j]<pop.getNumReal()) {
						double mass=pop.getMass(pool[n].parts[j]);
						pool[n].mass+=mass;
						pool[n].cmx+=mass*x[0];
						pool[n].cmy+=mass*x[1];
						pool[n].cmz+=mass*x[2];
					}
#endif
					for(int k=0; k<DIM; ++k) {
						if(x[k]<pool[n].bounds.min[k]) pool[n].bounds.min[k]=x[k];
						if(x[k]>pool[n].bounds.max[k]) pool[n].bounds.max[k]=x[k];
					}
				}
			} else {
				pool[n].maxRad=0.0;
				pool[n].velDispSum=0.0;
				pool[n].cnt=0.0;
			}
			if(pool[n].cnt>0) pool[n].searchRadius=5*sqrt(pool[n].velDispSum/pool[n].cnt)*pop.getTimeStep();
			else pool[n].searchRadius=0.0;
#ifdef GRAVITY
			if(pool[n].mass>0.0) {
				pool[n].cmx/=pool[n].mass;
				pool[n].cmy/=pool[n].mass;
				pool[n].cmz/=pool[n].mass;
				pool[n].size=pool[n].bounds.max[0]-pool[n].bounds.min[0];
				for(int i=1; i<DIM; ++i) {
					if(pool[n].size<pool[n].bounds.max[i]-pool[n].bounds.min[i])
						pool[n].size=pool[n].bounds.max[i]-pool[n].bounds.min[i];
				}
			}
#endif
		}

		template<class Population,class CollisionManager>
		void recurseForAll(Population &pop, CollisionManager &cm, int node1, KDBounds &b1, int node2, KDBounds &b2, std::vector<RecurseInfo> &info, int depth=0) {
			if (depth >= maxDepth)
			{
				info.resize(info.size()+1);
				info[info.size()-1].n1=node1;
				info[info.size()-1].n2=node2;
				info[info.size()-1].b1=b1;
				info[info.size()-1].b2=b2;
				return;
			}
			//printf("%d %d %e\n b1=(%e %e, %e %e)\n b2=(%e %e, %e %e)\n",node1,node2,searchRadius,b1.min[0],b1.max[0],b1.min[1],b1.max[1],b2.min[0],b2.max[0],b2.min[1],b2.max[1]);
			if(pool[node1].numParts>=0 && pool[node2].numParts>=0) {
				for(int i=0; i<pool[node1].numParts; ++i) {
					int j;
					if(node1==node2) j=i+1;
						else j=0;
					for(; j<pool[node2].numParts; ++j) {
						double t=pop.collisionTime(pool[node1].parts[i],pool[node2].parts[j]);
						if(t>=0.0 && t<pop.getTimeStep()) {
							cm.addPotentialCollision(pool[node1].parts[i],pool[node2].parts[j],t);
						}
					}
				}
			} else if(pool[node1].numParts>=0) {
				recurseForAll(pop,cm,node2,b2,node1,b1,info,depth+1);
			} else if(node1==node2) {
				// This is a special case that prevents me from
				// doing extra work. In this situation I know they
				// all overlap and I only need to recurse in 3 different
				// combination, not 4.
				int sd=pool[node1].splitDim;
				double sv=pool[node1].splitVal;
				double maxTmp=b1.max[sd];
				double minTmp=b1.min[sd];
				
				b1.max[sd]=b2.max[sd]=sv;
				recurseForAll(pop,cm,pool[node1].firstChild,b1,pool[node2].firstChild,b2,info,depth+1);
				b1.max[sd]=b2.max[sd]=maxTmp;
				b1.min[sd]=b2.min[sd]=sv;
				recurseForAll(pop,cm,pool[node1].firstChild+1,b1,pool[node2].firstChild+1,b2,info,depth+1);
				b1.min[sd]=minTmp;
				b1.max[sd]=sv;
				recurseForAll(pop,cm,pool[node1].firstChild,b1,pool[node2].firstChild+1,b2,info,depth+1);
				b1.max[sd]=maxTmp;
				b2.min[sd]=minTmp;
			} else {
				int sd=pool[node1].splitDim;
				double sv=pool[node1].splitVal;
				double tmp;
//				int pick=0;
				if(b2.min[sd]-(max(pool[node1].searchRadius,pool[node2].searchRadius)+pool[node1].maxRad+pool[node2].maxRad)<sv) {
					tmp=b1.max[sd];
					b1.max[sd]=sv;
					recurseForAll(pop,cm,node2,b2,pool[node1].firstChild,b1,info,depth+1);
					b1.max[sd]=tmp;
//					pick|=1;
				}
				if(b2.max[sd]+(max(pool[node1].searchRadius,pool[node2].searchRadius)+pool[node1].maxRad+pool[node2].maxRad)>sv) {
					tmp=b1.min[sd];
					b1.min[sd]=sv;
					recurseForAll(pop,cm,node2,b2,pool[node1].firstChild+1,b1,info,depth+1);
					b1.min[sd]=tmp;
//					pick|=2;
				}
//				printf("%d\n",pick);
			}
		}

		template<class Population,class CollisionManager>
		void recurseForOne(Population &pop,CollisionManager &cm,int p,int ignore,KDBounds &b,int node) {
			if(pool[node].numParts>=0) {
				for(int i=0; i<pool[node].numParts; ++i) {
					if(pool[node].parts[i]!=p && pool[node].parts[i]!=ignore) {
						double t=pop.collisionTime(p,pool[node].parts[i]);
						if(t>=0.0 && t<pop.getTimeStep()) {
							cm.addPotentialCollision(p,pool[node].parts[i],t);
						}
					}
				}
			} else {
				int sd=pool[node].splitDim;
				double sv=pool[node].splitVal;
				//double val=pop.get(p,sd);
				double val = location[p].c[sd];
				if(val-(pool[node].searchRadius+pool[node].maxRad+pop.getRadius(p))<sv) {
					double tmp=b.max[sd];
					b.max[sd]=sv;
					recurseForOne(pop,cm,p,ignore,b,pool[node].firstChild);
					b.max[sd]=tmp;
				}
				if(val+(pool[node].searchRadius+pool[node].maxRad+pop.getRadius(p))>sv){
					double tmp=b.min[sd];
					b.min[sd]=sv;
					recurseForOne(pop,cm,p,ignore,b,pool[node].firstChild+1);
					b.min[sd]=tmp;
				}
			}
		}
		
		template<class Population>
		void printTree(Population &pop,int node) {
			if(pool[node].numParts<0) {
				printf("Node %d %d %e\n",node,pool[node].splitDim,pool[node].splitVal);
				printf("left\n");
				printTree(pop,pool[node].firstChild);
				printf("right\n");
				printTree(pop,pool[node].firstChild+1);
				printf("done with %d\n",node);
			} else {
				printf("Particles in %d\n",node);
				for(int i=0; i<pool[node].numParts; ++i) {
					int p=pool[node].parts[i];
					printf("%d %e %e %e %e %e %e\n",p,pop.get(p,0),pop.get(p,1)
						,pop.get(p,2),pop.get(p,3),pop.get(p,4),pop.get(p,5));
				}
				printf("Done with %d\n",node);
			}
		}

		template<class Population>
		void doForce(int n,Population &pop,int i,double offsetX,double offsetY) {
			if(n<0) {
				printf("Negative node in force. %d %d\n",pc.getProcessNum(),n);
				exit(-1);
			}
			if(pool[n].mass<=0.0) return;
//			printf("%d %d %d\n",pc.getProcessNum(),n,i);
			if(pool[n].numParts>0) {
				for(int j=0; j<pool[n].numParts; ++j) {
					int oi=pool[n].parts[j];
					if(oi!=i && oi>=0 && oi<pop.getNumReal()) {
						double dx=pop.getx(oi)+offsetX-pop.getx(i);
						double dy=pop.gety(oi)+offsetY-pop.gety(i);
						double dz=pop.getz(oi)-pop.getz(i);
						double dist=sqrt(dx*dx+dy*dy+dz*dz);
						if(dist>(pop.getRadius(oi)+pop.getRadius(i))*0.9) {
							double mag=pop.getTimeStep()*pop.getMass(oi)/(dist*dist*dist);
							pop.setvx(i,pop.getvx(i)+dx*mag);
							pop.setvy(i,pop.getvy(i)+dy*mag);
							pop.setvz(i,pop.getvz(i)+dz*mag);
						}
					} else if(oi<0) {
//						printf("rp force %d of %d, %e %e %e %e\n",j,pool[n].numParts,remoteParticles[-oi].x,remoteParticles[-oi].y,remoteParticles[-oi].z,remoteParticles[-oi].mass);
						double dx=remoteParticles[-oi].x+offsetX-pop.getx(i);
						double dy=remoteParticles[-oi].y+offsetY-pop.gety(i);
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
			} else if(pool[n].numParts<0) {
				double dx=pool[n].cmx+offsetX-pop.getx(i);
				double dy=pool[n].cmy+offsetY-pop.gety(i);
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
					doForce(pool[n].firstChild,pop,i,offsetX,offsetY);
					doForce(pool[n].firstChild+1,pop,i,offsetX,offsetY);
				}
			}
		}

		/*
		 * This is a method that does gravity by doing a double recursion on two
		 * different nodes, the puller and the pulled.  The goal of this method is to
		 * reduce the amount of work spent doing distance calculations with top level
		 * nodes.
		 */
		template<class Population>
		void doPPForce(int pullerNode,int pulledNode,Population &pop,double offsetX,double offsetY,vector<GravRecurseInfo> &stack,int depth,vector<AccelVect> &acc) {
			if(depth>=maxDepth+3) {
				GravRecurseInfo gri(pullerNode,pulledNode,offsetX,offsetY);
				stack.push_back(gri);
				return;
			}
			if(pool[pullerNode].mass<=0.0 || pool[pulledNode].mass<=0.0) return;
			if(pool[pullerNode].numParts>0 || pool[pulledNode].numParts>0) {
				if(depth>=0) {
					GravRecurseInfo gri(pullerNode,pulledNode,offsetX,offsetY);
					stack.push_back(gri);
					return;
				}
				pullAllChildren(pullerNode,pulledNode,pop,offsetX,offsetY,acc);
			} else {
				if(pullerNode==pulledNode && offsetX==0.0 && offsetY==0) {
					doPPForce(pool[pullerNode].firstChild,pool[pullerNode].firstChild,pop,offsetX,offsetY,stack,depth+1,acc);
					doPPForce(pool[pullerNode].firstChild,pool[pullerNode].firstChild+1,pop,offsetX,offsetY,stack,depth+2,acc);
					doPPForce(pool[pullerNode].firstChild+1,pool[pullerNode].firstChild,pop,offsetX,offsetY,stack,depth+2,acc);
					doPPForce(pool[pullerNode].firstChild+1,pool[pullerNode].firstChild+1,pop,offsetX,offsetY,stack,depth+1,acc);
				} else {
					double ndx=pool[pullerNode].cmx+offsetX-pool[pulledNode].cmx;
					double ndy=pool[pullerNode].cmy+offsetY-pool[pulledNode].cmy;
					double ndz=pool[pullerNode].cmz-pool[pulledNode].cmz;
					double ndist=sqrt(ndx*ndx+ndy*ndy+ndz*ndz);
					if(theta*(ndist-pool[pulledNode].size)>pool[pullerNode].size) {
						// Far away so pull all
						if(depth>=0) {
							GravRecurseInfo gri(pullerNode,pulledNode,offsetX,offsetY);
							stack.push_back(gri);
							return;
						}
						pullAllChildren(pullerNode,pulledNode,pop,offsetX,offsetY,acc);
					} else if(theta*(ndist+pool[pulledNode].size)<pool[pullerNode].size) {
						// Nearby so recurse on the puller
						doPPForce(pool[pullerNode].firstChild,pulledNode,pop,offsetX,offsetY,stack,depth+2,acc);
						doPPForce(pool[pullerNode].firstChild+1,pulledNode,pop,offsetX,offsetY,stack,depth+2,acc);
					} else {
						// in between so recurse on the pulled
						doPPForce(pullerNode,pool[pulledNode].firstChild,pop,offsetX,offsetY,stack,depth+3,acc);
						doPPForce(pullerNode,pool[pulledNode].firstChild+1,pop,offsetX,offsetY,stack,depth+3,acc);
					}
				}
			}
		}

		/*
		 * Method to run down the pulledNode to all leaves and apply gravity from the pullerNode.
		 */
		template<class Population>
		void pullAllChildren(int pullerNode,int pulledNode,Population &pop,double offsetX,double offsetY,vector<AccelVect> &acc) {
			if(pool[pulledNode].numParts>0) {
				if(pool[pullerNode].numParts>0) {
					int numk=pool[pullerNode].numParts;
					for(int k=0; k<numk; ++k) {
						int oi=pool[pullerNode].parts[k];
						if(oi>=0 && oi<pop.getNumReal()) {
							int j=0;
							if(pullerNode==pulledNode) j=k+1;
							for(; j<pool[pulledNode].numParts; ++j) {
								int i=pool[pulledNode].parts[j];
								double dx=pop.getx(oi)+offsetX-pop.getx(i);
								double dy=pop.gety(oi)+offsetY-pop.gety(i);
								double dz=pop.getz(oi)-pop.getz(i);
								double dist=sqrt(dx*dx+dy*dy+dz*dz);
								if(dist==0.0) {
									printf("ERROR Particle: %d %d %d %d %d\n",pullerNode,pulledNode,i,oi,j);
								}
								double mag=pop.getTimeStep()*pop.getMass(oi)/(dist*dist*dist);
								acc[i].ax+=dx*mag;
								acc[i].ay+=dy*mag;
								acc[i].az+=dz*mag;
								if(pullerNode==pulledNode) {
									mag=-pop.getTimeStep()*pop.getMass(i)/(dist*dist*dist);
									acc[oi].ax+=dx*mag;
									acc[oi].ay+=dy*mag;
									acc[oi].az+=dz*mag;
								}
							}
						} else if(oi<0) {
							for(int j=0; j<pool[pulledNode].numParts; ++j) {
								int i=pool[pulledNode].parts[j];
								double dx=remoteParticles[-oi].x+offsetX-pop.getx(i);
								double dy=remoteParticles[-oi].y+offsetY-pop.gety(i);
								double dz=remoteParticles[-oi].z-pop.getz(i);
								double dist=sqrt(dx*dx+dy*dy+dz*dz);
								if(dist==0.0) {
									printf("ERROR Node: %d %d %d %d %d\n",pullerNode,pulledNode,i,oi,j);
								}
								if(dist>2*pop.getRadius(i)) {
									double mag=pop.getTimeStep()*remoteParticles[-oi].mass/(dist*dist*dist);
									acc[i].ax+=dx*mag;
									acc[i].ay+=dy*mag;
									acc[i].az+=dz*mag;
								}
							}
						}
					}
				} else {
					for(int k=0; k<pool[pulledNode].numParts; ++k) {
						int i=pool[pulledNode].parts[k];
						double dx=pool[pullerNode].cmx+offsetX-pop.getx(i);
						double dy=pool[pullerNode].cmy+offsetY-pop.gety(i);
						double dz=pool[pullerNode].cmz-pop.getz(i);
						double dsqr=dx*dx+dy*dy+dz*dz;
						double dist=sqrt(dsqr);
						if(dist==0.0) {
							printf("ERROR Node 2: %d %d %d %d\n",pullerNode,pulledNode,i,k);
						}
						double mag=pop.getTimeStep()*pool[pullerNode].mass/(dist*dsqr);
						acc[i].ax+=dx*mag;
						acc[i].ay+=dy*mag;
						acc[i].az+=dz*mag;
					}
				}
			} else if(pool[pulledNode].numParts<0){
				pullAllChildren(pullerNode,pool[pulledNode].firstChild,pop,offsetX,offsetY,acc);
				pullAllChildren(pullerNode,pool[pulledNode].firstChild+1,pop,offsetX,offsetY,acc);
			}
		}

		template<class Population>
		void communicateTrees(Population &pop) {
			receivedTrees.resize(0);
			remoteParticles.resize(1);
			vector<int> rt;
			for(int i=0; i<pc.getNumProcesses(); ++i) {
				int numRead;
				if(i==pc.getProcessNum()) {
					// receive from others.
					for(int j=0; j<DIM; ++j) {
						boundsBuffer[2*j]=pool[0].bounds.min[j];
						boundsBuffer[2*j+1]=pool[0].bounds.max[j];
					}
					for(int j=0; j<pc.getNumProcesses(); ++j) {
						if(j!=i) {
							//printf("%d recieve from %d\n",i,j);
							pc.sendTo(j,boundsBuffer);
							pc.readFrom(j,sizeBuffer,numRead);
							treeBuffer.resize((int)sizeBuffer[0]);
							pc.readFrom(j,treeBuffer,numRead);
							rt.resize(0);
							//printf("%d reconstruct\n",i);
							reconstructTree(rt);
							for(unsigned int k=0; k<rt.size(); ++k) {
								receivedTrees.push_back(rt[k]);
							}
							//printf("%d going to next\n",i);
						}
					}
				} else {
					// send to others.
//					printf("%d send to %d\n",pc.getProcessNum(),i);
					pc.readFrom(i,boundsBuffer,numRead);
//					printf("%d pack tree\n",pc.getProcessNum());
					treeBuffer.resize(0);
					packTree(pop,boundsBuffer[0],boundsBuffer[1],boundsBuffer[2],boundsBuffer[3],boundsBuffer[4],boundsBuffer[5]);
//					printf("%d done packing tree\n",pc.getProcessNum());
#ifdef AZIMUTHAL_MIRRORS
					double cellWidth=bounds.getMaxYTotal()-bounds.getMinYTotal();
					packTree(pop,boundsBuffer[0],boundsBuffer[1],boundsBuffer[2],boundsBuffer[3],boundsBuffer[4],boundsBuffer[5],0,cellWidth);
					packTree(pop,boundsBuffer[0],boundsBuffer[1],boundsBuffer[2],boundsBuffer[3],boundsBuffer[4],boundsBuffer[5],0,-cellWidth);
#endif
#ifdef FULL_MIRRORS
		            double shearOffset=bounds.getShearOffset();
		            double cellHeight=bounds.getMaxX()-bounds.getMinX();
					packTree(pop,boundsBuffer[0],boundsBuffer[1],boundsBuffer[2],boundsBuffer[3],boundsBuffer[4],boundsBuffer[5],cellHeight,shearOffset+cellWidth);
	                packTree(pop,boundsBuffer[0],boundsBuffer[1],boundsBuffer[2],boundsBuffer[3],boundsBuffer[4],boundsBuffer[5],cellHeight,shearOffset);
	                packTree(pop,boundsBuffer[0],boundsBuffer[1],boundsBuffer[2],boundsBuffer[3],boundsBuffer[4],boundsBuffer[5],cellHeight,shearOffset-cellWidth);
	                packTree(pop,boundsBuffer[0],boundsBuffer[1],boundsBuffer[2],boundsBuffer[3],boundsBuffer[4],boundsBuffer[5],-cellHeight,-shearOffset+cellWidth);
	                packTree(pop,boundsBuffer[0],boundsBuffer[1],boundsBuffer[2],boundsBuffer[3],boundsBuffer[4],boundsBuffer[5],-cellHeight,-shearOffset);
	                packTree(pop,boundsBuffer[0],boundsBuffer[1],boundsBuffer[2],boundsBuffer[3],boundsBuffer[4],boundsBuffer[5],-cellHeight,-shearOffset-cellWidth);
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
		 * 9 - minz
		 * 10 - maxz
		 * 11 - splitDim
		 * 12 - spitVal
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
		void packTree(Population &pop,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double offsetX=0.0,double offsetY=0.0) {
			packNode(pop,offsetX,offsetY,0);
			packRecur(pop,xmin,xmax,ymin,ymax,zmin,zmax,offsetX,offsetY,0);
		}
		
		template<class Population>
		void packRecur(Population &pop,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double offsetX,double offsetY,int n) {
//			printf("%d pack recur on %d\n",pc.getProcessNum(),n);
			if(pool[n].numParts<1) {
				// Now check distance and see if we recurse.
				double xSep=std::max(xmin-(pool[n].bounds.max[0]+offsetX),(pool[n].bounds.min[0]+offsetX)-xmax);
				if(xSep<0.0) xSep=0.0;
				double ySep=std::max(ymin-(pool[n].bounds.max[1]+offsetY),(pool[n].bounds.min[1]+offsetY)-ymax);
				if(ySep<0.0) ySep=0.0;
				double zSep=std::max(zmin-(pool[n].bounds.max[2]),(pool[n].bounds.min[2])-zmax);
				if(zSep<0.0) zSep=0.0;
				if((xSep*xSep+ySep*ySep+zSep*zSep)*theta*theta<=pool[n].size*pool[n].size) {
					packNode(pop,offsetX,offsetY,pool[n].firstChild);
					packNode(pop,offsetX,offsetY,pool[n].firstChild+1);
					packRecur(pop,xmin,xmax,ymin,ymax,zmin,zmax,offsetX,offsetY,pool[n].firstChild);
					packRecur(pop,xmin,xmax,ymin,ymax,zmin,zmax,offsetX,offsetY,pool[n].firstChild+1);
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
				treeBuffer.push_back(pool[n].bounds.min[0]+offsetX);
				treeBuffer.push_back(pool[n].bounds.max[0]+offsetX);
				treeBuffer.push_back(pool[n].bounds.min[1]+offsetY);
				treeBuffer.push_back(pool[n].bounds.max[1]+offsetY);
				treeBuffer.push_back(pool[n].bounds.min[2]);
				treeBuffer.push_back(pool[n].bounds.max[2]);
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
				//printf("%d reconstruct %d of %d\n",pc.getProcessNum(),i,treeBuffer.size());
				if(firstFree>=pool.size()) {
					pool.resize(pool.size()+pool.size()/2);
				}
				pool[firstFree].firstChild=-1;
				if(treeBuffer[i]==0) {
					//printf("CM node\n");
					pool[firstFree].cmx=treeBuffer[i+1];
					pool[firstFree].cmy=treeBuffer[i+2];
					pool[firstFree].cmz=treeBuffer[i+3];
					pool[firstFree].mass=treeBuffer[i+4];
					pool[firstFree].bounds.min[0]=treeBuffer[i+5];
					pool[firstFree].bounds.max[0]=treeBuffer[i+6];
					pool[firstFree].bounds.min[1]=treeBuffer[i+7];
					pool[firstFree].bounds.max[1]=treeBuffer[i+8];
					pool[firstFree].bounds.min[2]=treeBuffer[i+9];
					pool[firstFree].bounds.max[2]=treeBuffer[i+10];
					pool[firstFree].splitDim=(int)treeBuffer[i+11];
					pool[firstFree].splitVal=treeBuffer[i+12];
					pool[firstFree].numParts=-1;
					i+=13;
				} else if(treeBuffer[i]==1) {
					pool[firstFree].numParts=(int)treeBuffer[i+1];
					//printf("part node %d\n",pool[firstFree].numParts);
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
							pool[firstFree].bounds.min[0]=x;
							pool[firstFree].bounds.max[0]=x;
							pool[firstFree].bounds.min[1]=y;
							pool[firstFree].bounds.max[1]=y;
							pool[firstFree].bounds.min[2]=z;
							pool[firstFree].bounds.max[2]=z;
						} else {
							if(x<pool[firstFree].bounds.min[0]) pool[firstFree].bounds.min[0]=x;
							if(x>pool[firstFree].bounds.max[0]) pool[firstFree].bounds.max[0]=x;
							if(y<pool[firstFree].bounds.min[1]) pool[firstFree].bounds.min[1]=y;
							if(y>pool[firstFree].bounds.max[1]) pool[firstFree].bounds.max[1]=y;
							if(z<pool[firstFree].bounds.min[2]) pool[firstFree].bounds.min[2]=z;
							if(z>pool[firstFree].bounds.max[2]) pool[firstFree].bounds.max[2]=z;
						}
					}
				} else {
					printf("Error: Bad code in treeBuffer! %d is %f\n",i,treeBuffer[i]);
					exit(-1);
				}
//				printf("Bounds %e %e %e %e\n",
//					pool[firstFree].minx,pool[firstFree].maxx,
//					pool[firstFree].miny,pool[firstFree].maxy);
				pool[firstFree].size=std::max(std::max(pool[firstFree].bounds.max[0]-pool[firstFree].bounds.min[0],pool[firstFree].bounds.max[1]-pool[firstFree].bounds.min[1]),pool[firstFree].bounds.max[2]-pool[firstFree].bounds.min[2]);
				if(firstFree!=curRoot) {
					if(fitsIn(firstFree,curRoot)) {
						if(addCnt==1) {
							//printf("Add node %d %d\n",curRoot,firstFree);
							addNode(curRoot,firstFree);
							addCnt=0;
							//printf("Done adding\n");
						} else {
							addCnt=1;
						}
					} else {
						if(addCnt==0) {
							printf("Error: new root when expecting second child.\n");
							printf("new bounds %e %e %e %e %e %e\n",pool[firstFree].bounds.min[0],pool[firstFree].bounds.max[0],pool[firstFree].bounds.min[1],pool[firstFree].bounds.max[1],pool[firstFree].bounds.min[2],pool[firstFree].bounds.max[2]);
							printf("root bounds %e %e %e %e %e %e\n",pool[curRoot].bounds.min[0],pool[curRoot].bounds.max[0],pool[curRoot].bounds.min[1],pool[curRoot].bounds.max[1],pool[curRoot].bounds.min[2],pool[curRoot].bounds.max[2]);
							exit(-1);
						}
//						printf("%d Diff tree %d %d (%e %e) (%e %e) (%e %e) (%e %e)\n",pc.getProcessNum(),curRoot,firstFree,
//							pool[firstFree].minx,pool[curRoot].minx,
//							pool[firstFree].maxx,pool[curRoot].maxx,
//							pool[firstFree].miny,pool[curRoot].miny,
//							pool[firstFree].maxy,pool[curRoot].maxyi);
						curRoot=firstFree;
						roots.push_back(curRoot);
						addCnt=1;
					}
				} else {
//					printf("%d Diff tree %d %e %e %e %e\n",pc.getProcessNum(),firstFree,
//						pool[firstFree].minx,
//						pool[firstFree].maxx,
//						pool[firstFree].miny,
//						pool[firstFree].maxy);
				}
				firstFree++;
			}
		}
		
		/**
		 * Checks if node n fits inside the bounds of node m.
		 */
		bool fitsIn(int n,int m) {
			return pool[n].bounds.min[0]>=pool[m].bounds.min[0] &&
				pool[n].bounds.max[0]<=pool[m].bounds.max[0] &&
				pool[n].bounds.min[1]>=pool[m].bounds.min[1] &&
				pool[n].bounds.max[1]<=pool[m].bounds.max[1] &&
				pool[n].bounds.min[2]>=pool[m].bounds.min[2] &&
				pool[n].bounds.max[2]<=pool[m].bounds.max[2];
		}
		
		void addNode(int n,int newNode) {
			int p=-1;
			double cx=(pool[newNode].bounds.max[0]+pool[newNode].bounds.min[0])*0.5;
			double cy=(pool[newNode].bounds.max[1]+pool[newNode].bounds.min[1])*0.5;
			double cz=(pool[newNode].bounds.max[2]+pool[newNode].bounds.min[2])*0.5;
			while(n!=-1) {
				//printf("%d %d %d\n",pc.getProcessNum(),n,p);
				p=n;
				n=pool[n].childNode(cx,cy,cz);
				if(n==0) {
					printf("Error %d: add node got to 0, p=%d, newNode=%d, pool.size=%d.\n",pc.getProcessNum(),p,newNode,pool.size());
					printf("p.numParts=%d, n.numParts=%d\n",pool[p].numParts,pool[newNode].numParts);
					printf("cx=%e, cy=%e\n",cx,cy);
					printf("p.splitDim=%d, p.splitVal=%e\n",pool[p].splitDim,pool[p].splitVal);
					printf("p bounds %e %e %e %e %e %e\n",pool[p].bounds.min[0],pool[p].bounds.max[0],pool[p].bounds.min[1],pool[p].bounds.max[1],pool[p].bounds.min[2],pool[p].bounds.max[2]);
					printf("new bounds %e %e %e %e %e %e\n",pool[newNode].bounds.min[0],pool[newNode].bounds.max[0],pool[newNode].bounds.min[1],pool[newNode].bounds.max[1],pool[newNode].bounds.min[2],pool[newNode].bounds.max[2]);
					if(pool[newNode].numParts>0) {
						for(int i=0; i<pool[newNode].numParts; ++i) {
							int part=pool[newNode].parts[i];
							printf("Part %d at %e %e %e\n",part,remoteParticles[-part].x,remoteParticles[-part].y,remoteParticles[-part].z);
						}
					}
					fflush(stdout);
					exit(-1);
				}
			}
			pool[p].firstChild=newNode;
			//printf("Putting %d under %d.\n",newNode,p);
			if(!fitsIn(newNode,p)) {
				printf("Put %d in %d but it doesn't fit.\n",newNode,p);
				printf("p bounds %e %e %e %e\n",pool[p].bounds.min[0],pool[p].bounds.max[0],pool[p].bounds.min[1],pool[p].bounds.max[1]);
				printf("new bounds %e %e %e %e\n",pool[newNode].bounds.min[0],pool[newNode].bounds.max[0],pool[newNode].bounds.min[1],pool[newNode].bounds.max[1]);
			}
		}
		
		std::vector<KDNode> pool;
		std::vector<Location> location;
		std::vector<std::vector<bool> > inUse;
		vector<vector<AccelVect> > acc;
		double gridSpacing;
		int maxDepth;
		int firstFree;
		KDBounds totalBounds;
		double theta;
		ProcessorCommunication &pc;
        Boundary &bounds;  // boundary conditions
		double minDist;
		omp_lock_t lock;
		vector<double> boundsBuffer;
		vector<double> sizeBuffer;
		vector<double> treeBuffer;
		vector<int> receivedTrees;
		vector<RemoteParticle> remoteParticles;
};

#endif
