// This is my version of a spatial KD tree used for collision detection and
// gravity.  Note that if this is used for gravity and collisions in a
// DoubleForce, the gravity needs to come first.

#ifndef GRAV_COLL_TREE
#define GRAV_COLL_TREE

#include <vector>
#include <algorithm>
#include <unordered_set>
#include <tuple>
#include <functional>

#include "ParticleIndex.h"
#include "ShortRangeForces.h"
#include "AccelVect.h"

//#define LOCK_GRAPH

const int MAX_NUM=2;

const int DIM=3;

using std::vector;
using std::min;
using std::max;
using std::unordered_set;
using std::tuple;
using std::make_tuple;

/* Uncomment after fixing ParticleIndex
struct KDNodeIndex {
	int i;
	KDNodeIndex operator+(int diff) {
		return KDNodeIndex{i+diff};
	}
}
*/
#ifdef GRAV_STATS
class IterCounter {
	public:
		static void recordSqr() {
			sqrCnt++;
		}
		static void recordSqrt() {
			sqrtCnt++;
		}
		static void printCounts(FILE *fout) {
			fprintf(fout,"sqr=%ld\n",sqrCnt);
			fprintf(fout,"sqrt=%ld\n",sqrtCnt);
		}

		static long long int sqrCnt;
		static long long int sqrtCnt;
};

long long int IterCounter::sqrCnt;
long long int IterCounter::sqrtCnt;

#endif

struct TwoTupleHash {
	size_t operator()(const tuple<int, int> &t) const {
		return std::get<0>(t) + 13*std::get<1>(t);
	}
};

struct CollisionSetBuilder {
	unordered_set<tuple<int, int>, TwoTupleHash> &treeSet;
	CollisionSetBuilder(unordered_set<tuple<int, int>, TwoTupleHash> &ts): treeSet(ts) {}

	template<typename Population>
	void operator()(Population &pop, ParticleIndex pi, ParticleIndex pj) {
		double dx = pop.getx(pi) - pop.getx(pj);
		double dy = pop.gety(pi) - pop.gety(pj);
		double dz = pop.getz(pi) - pop.getz(pj);
		if (dx*dx + dy*dy + dz*dz < (pop.getRadius(pi)+pop.getRadius(pj)) * (pop.getRadius(pi)+pop.getRadius(pj))) {
			int pmin = min(pi.i, pj.i);
			int pmax = max(pi.i, pj.i);
			treeSet.insert(make_tuple(pmin, pmax));
		}
	}
};

struct KDNode {
	int splitDim;
	double splitVal;
	int firstChild;
	ParticleIndex parts[MAX_NUM];
	int numParts;
	
	double cmx,cmy,cmz,mass;
	double min[DIM],max[DIM],size;
	double maxRad;
	double cnt;
	double velDispSum;
	double searchRadius;
#ifdef LOCK_GRAPH
	bool inUse;
#endif
	
	int childNode(double x,double y,double z) {
		if((splitDim==0 && x<=splitVal) || (splitDim==1 && y<=splitVal) || (splitDim==2 && z<=splitVal)) return firstChild;
		return firstChild+1;
	}

	bool operator<(const KDNode &n) { return min[1]<n.min[1]; }
};

struct Location {
	double c[DIM];
};

#ifdef PARALLEL
struct RemoteParticle {
	RemoteParticle() {}
	RemoteParticle(double nx,double ny,double nz,double nmass):x(nx),y(ny),z(nz),mass(nmass) {}
	double x,y,z,mass;
};
#endif

struct RecurseInfo {
	int n1,n2;
	int end;
};

struct GravRecurseInfo {
	GravRecurseInfo(int n1,int n2,double ox,double oy):puller(n1),pulled(n2),offsetX(ox),offsetY(oy) {}
	int puller,pulled;
	double offsetX,offsetY;
};

template<class Boundary, class ShortRangeForce = ShortRangeGravityOnly>
class GravCollTree {
	public:
#ifndef PARALLEL
		GravCollTree(double t,Boundary &b, ShortRangeForce srf = ShortRangeGravityOnly(),double md=0):theta(t),bounds(b),shortRangeForce(srf),minDist(md) {
#ifdef PFORCE
			acc.resize(omp_get_max_threads());
#else
			acc.resize(1);
#endif
			maxDepth = (int) (log(omp_get_max_threads())/log(2) + 2);
			omp_init_lock(&lock);
#ifdef REORDER
			stepCount=0;
#endif
		}
#else
		GravCollTree(double t,ProcessorCommunication &procComm,Boundary &b, ShortRangeForce srf = ShortRangeGravityOnly(),double md=0):theta(t),pc(procComm),bounds(b),shortRangeForce(srf),minDist(md),boundsBuffer(6),sizeBuffer(1) {
#ifdef PFORCE
			acc.resize(omp_get_max_threads());
#else
			acc.resize(1);
#endif
			maxDepth = (int) (log(omp_get_max_threads())/log(2) + 2);
			omp_init_lock(&lock);
		}
#endif

		~GravCollTree() {
			omp_destroy_lock(&lock);
		}

		template<class Population>
		void applyForce(Population &pop) {
			const int nt=acc.size();
			const int nb=pop.getNumBodies();
			#pragma omp parallel for
			for(int i=0; i<nt; ++i) {
				acc[i].resize(pop.getNumBodies());
				for(int j=0; j<nb; ++j) {
					acc[i][j].ax=0.0;
					acc[i][j].ay=0.0;
					acc[i][j].az=0.0;
				}
			}
			//printf("Do tree gravity.\n");
			buildReal(pop);
#ifdef REORDER
			if(stepCount%10==0) {
				vector<int> map(pop.getNumBodies());
				vector<int> invmap(pop.getNumBodies());
				for(int i=0; i<map.size(); ++i) {
					map[i]=i;
					invmap[i]=i;
				}
				int nextPos=0;
				reorderFromTree(pop,0,map,invmap,nextPos);
			}
			stepCount++;
#endif
#ifdef PARALLEL
			const double cellWidth=bounds.getMaxYTotal()-bounds.getMinYTotal();
#else
			const double cellWidth=bounds.getMaxY()-bounds.getMinY();
#endif
#ifdef FULL_MIRRORS
            const double shearOffset=bounds.getShearOffset();
            const double cellHeight=bounds.getMaxX()-bounds.getMinX();
#endif
#ifdef PARALLEL
			communicateTrees(pop);
#endif
#ifdef PFORCE
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
#ifdef PARALLEL
			for(unsigned int j=0; j<receivedTrees.size(); ++j) {
				doPPForce(receivedTrees[j],0,pop,0,0,stack,0,acc[0]);
			}
#endif
			int ssize = (int) stack.size();
			#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < ssize; ++i) {
				doPPForce(stack[i].puller,stack[i].pulled,pop,stack[i].offsetX,stack[i].offsetY,stack,-1000000,acc[omp_get_thread_num()]);
			}
#else // for PFORCE
#ifdef DEBUG
			verifyOverlaps(pop);
#endif
			#pragma omp parallel for schedule(dynamic,100)
			for(int ii=0; ii<nb; ++ii) {
				ParticleIndex pi = {ii};
				doForce(0, pop, pi, 0, 0, acc[0]);
#ifdef AZIMUTHAL_MIRRORS
				for(int j=1; j<=AZIMUTHAL_MIRRORS; ++j) {
	                doForce(0,pop,pi,0,cellWidth*j, acc[0]);
                	doForce(0,pop,pi,0,-cellWidth*j, acc[0]);
#ifdef FULL_MIRRORS
					for(int k=1; k<=FULL_MIRRORS; ++k) {
						double totalShear=k*shearOffset;
						while(totalShear<bounds.getMinY()) totalShear+=cellWidth;
	                	doForce(0,pop,pi,cellHeight*k,totalShear+cellWidth*j, acc[0]);
	                	doForce(0,pop,pi,cellHeight*k,totalShear-cellWidth*j, acc[0]);
	                	doForce(0,pop,pi,-cellHeight*k,-totalShear+cellWidth*j, acc[0]);
	                	doForce(0,pop,pi,-cellHeight*k,-totalShear-cellWidth*j, acc[0]);
					}
#endif
				}
#ifdef FULL_MIRRORS
				for(int k=1; k<=FULL_MIRRORS; ++k) {
					double totalShear=k*shearOffset;
					while(totalShear<bounds.getMinY()) totalShear+=cellWidth;
                	doForce(0,pop,pi,cellHeight*k,totalShear, acc[0]);
                	doForce(0,pop,pi,-cellHeight*k,-totalShear, acc[0]);
				}
#endif
#endif
#ifdef PARALLEL
				for(unsigned int j=0; j<receivedTrees.size(); ++j) {
					doForce(receivedTrees[j],pop,pi,0,0, acc[0]);
				}
#endif
			}
#endif
			#pragma omp parallel for
			for(int jj=0; jj<nb; ++jj) {
				ParticleIndex pj = {jj};
				double ax = 0.0;
				double ay = 0.0;
				double az = 0.0;
				for(unsigned int i=0; i<acc.size(); ++i) {
					ax += acc[i][pj.i].ax;
					ay += acc[i][pj.i].ay;
					az += acc[i][pj.i].az;
				}
				if (jj == 0) printf("Acc = %e %e %e\n", ax, ay, az);
				pop.setvx(pj,pop.getvx(pj)+pop.getTimeStep() * ax);
				pop.setvy(pj,pop.getvy(pj)+pop.getTimeStep() * ay);
				pop.setvz(pj,pop.getvz(pj)+pop.getTimeStep() * az);
				pop.adjustAfterForce(pj);
			}
		}

		template<class Population>
		void build(Population &pop) {
#ifndef GRAVITY
			buildReal(pop);
#endif
		}

		template<class Population>
		void buildReal(Population &pop) {
			//printf("Building.\n");
			pool.resize(pop.getNumBodies()*2);
#ifdef LOCK_GRAPH
			int psize=pool.size();
			#pragma omp parallel for
			for(int i=0; i<psize; ++i) {
				pool[i].inUse=false;
			}
#endif
			location.resize(pop.getNumBodies());
			pool[0].parts[0]={0};
			pool[0].numParts=0;
			firstFree=1;
			const int nb=pop.getNumBodies();
			#pragma omp parallel for
			for(int ii=0; ii<nb; ++ii) {
				ParticleIndex pi = {ii};
				for(int j=0; j<DIM; ++j) {
					location[pi.i].c[j] = pop.get(pi,j);
				}
			}

			vector<ParticleIndex> index(pop.getNumBodies());
			#pragma omp parallel for
			for(int i=0; i<nb; ++i) {
				index[i]=ParticleIndex{i};
			}

			//printf("Build recur to %d.\n",maxDepth);
			fflush(stdout);
			vector<RecurseInfo> info;
			buildRecur(pop,index,0,pop.getNumReal(),0,info);
			int isize = (int) info.size();
			//printf("Build loop. %d\n",isize);
			fflush(stdout);
			#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < isize; i++)
				buildRecur(pop,index,info[i].n2,info[i].end,info[i].n1,info,-1000000);

			info.resize(0);
			//printf("finalize top.\n");
			fflush(stdout);
			finalize(0,pop,info);
			isize = (int) info.size();
			//printf("finalize loop. %d\n",isize);
			fflush(stdout);
			#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < isize; i++)
				finalize(info[i].n1,pop,info,-1000000);
			info.resize(0);
			//printf("finalize end.\n");
			fflush(stdout);
			finalize(0,pop,info);

			printf("Search radius at top is %e\n",2*pool[0].maxRad+pool[0].searchRadius);

			//printf("End build %d nodes\n",firstFree);
			fflush(stderr);
//			printTree(pop,root);
		}

		double getMinSpacing() {
			return 3*(2*pool[0].maxRad+pool[0].searchRadius);
		}

		template<class Population,class ActionType>
		void findAllCollisions(Population &pop,ActionType &action) {
#ifndef LOCK_GRAPH
			buildLockGrid(pop);
#else
			buildLockGraph(pop);
#endif

			runThroughCollisionPairs(pop,action);
		}

		template<class Population,typename ActionType>
		void runThroughCollisionPairs(Population &pop, ActionType& action) {
			//Recurse through the top maxDepth levels of the tree.
			//Instead of making further recursive calls, store the parameters for
			//each call as a RecurseInfo and store it in the info queue.
			std::vector<RecurseInfo> info(0);
			recurseForAll(pop,action,0,0,info);
			
			//Make all stored recursive calls in parallel.
			int isize = (int) info.size();
			#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < isize; i++)
				recurseForAll(pop,action,info[i].n1,info[i].n2,info,-1000000);
		}

		template<class Population,class CollisionManager>
		void findCollisions(Population &pop,CollisionManager &cm,ParticleIndex p,ParticleIndex ignore) {
			recurseForOne(pop,cm,p,ignore,0);
		}

#ifndef LOCK_GRAPH
		int getBinX(ParticleIndex p) {
			return (int) ((location[p.i].c[0] - pool[0].min[0])/gridSpacing);
		}

		int getBinY(ParticleIndex p) {
			return (int) ((location[p.i].c[1] - pool[0].min[1])/gridSpacing);
		}

		bool isSafe(ParticleIndex pi) {
			int x=getBinX(pi);
			int y=getBinY(pi);
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

		void setInUse(ParticleIndex pi,bool val) {
			int x=getBinX(pi);
			int y=getBinY(pi);
			inUse[x][y] = val;
		}
#else
		bool isSafe(ParticleIndex pi) {
			int lockList=particleLockNodes[pi.i];
			if(pool[lockList].inUse) return false;
			int num=nodeEdges[lockList].size();
			for(int j=0; j<num; ++j) {
				if(pool[nodeEdges[lockList][j]].inUse) return false;
			}
			return true;
		}

		void setInUse(ParticleIndex pi,bool val) {
			pool[particleLockNodes[pi.i]].inUse = val;
		}
#endif

		template<class Population>
		bool verifyOverlaps(Population &pop) {
			unordered_set<tuple<int, int>, TwoTupleHash> bfSet;
			unordered_set<tuple<int, int>, TwoTupleHash> treeSet;

			for (ParticleIndex pi = {0}; pi < pop.getNumBodies(); ++pi.i) {
				for (ParticleIndex pj = {pi.i+1}; pj < pop.getNumBodies(); ++pj.i) {
					double dx = pop.getx(pi) - pop.getx(pj);
					double dy = pop.gety(pi) - pop.gety(pj);
					double dz = pop.getz(pi) - pop.getz(pj);
					if (dx*dx + dy*dy + dz*dz < (pop.getRadius(pi)+pop.getRadius(pj)) * (pop.getRadius(pi)+pop.getRadius(pj))) {
						bfSet.insert(make_tuple(pi.i, pj.i));
					}
				}
			}
			CollisionSetBuilder action(treeSet);
			runThroughCollisionPairs(pop, action);

			printf("BF overlaps: %ld, Tree overlaps: %ld\n", bfSet.size(), treeSet.size());
			if (bfSet != treeSet) {
				for (auto pair: treeSet) bfSet.erase(pair);
				printf("Missing overlaps: %ld\n", bfSet.size());
				for (auto pair: bfSet) printf("%d %d\n", std::get<0>(pair), std::get<1>(pair));
			}

			return bfSet == treeSet;
		}



	private:
#ifndef LOCK_GRAPH
		//To synchronize the KDTree we must create a grid back-end to keep track
		//of which regions contain particles that are currently being processed.
		template<class Population>
		void buildLockGrid(Population &pop) {
			gridSpacing = 2*pool[0].maxRad + pool[0].searchRadius;
			int xlen = (int) ((pool[0].max[0] - pool[0].min[0])/gridSpacing + 1);
			int ylen = (int) ((pool[0].max[1] - pool[0].min[1])/gridSpacing + 1);
			printf("%d (%e %e %e) %d (%e %e)\n",xlen,pool[0].max[0],pool[0].min[0],gridSpacing,ylen,pool[0].max[1],pool[0].min[1]);
			inUse.resize(xlen);
			for (int i = 0; i < xlen; i++)
			{
				inUse[i].resize(ylen);
				for (int j = 0; j < ylen; j++)
					inUse[i][j] = false;
			}
		}
#else
		template<class Population>
		void buildLockGraph(Population &pop) {
			particleLockNodes.resize(pop.getNumReal());
			nodeEdges.resize(pool.size());
			for(unsigned int i=0; i<nodeEdges.size(); ++i) {
				nodeEdges[i].resize(0);
				pool[i].inUse=false;
			}
			vector<RecurseInfo> info(0);
			buildLockGraphRecur(0,0,info);
			int isize = (int) info.size();
			#pragma omp parallel for schedule(dynamic)
			for (int i=0; i<isize; ++i)
				buildLockGraphRecur(info[i].n1,info[i].n2,info,-1000000);
//			int cnt=0,sum=0,leafCnt=0;
//			for(int i=0; i<nodeEdges.size(); ++i) {
//				if(nodeEdges[i].size()>0) cnt++;
//				if(pool[i].numParts>0) leafCnt++;
//				sum+=nodeEdges[i].size();
//			}
//			printf("%d lock nodes, %d leaf nodes\n",cnt,leafCnt);
//			printf("Average links=%f\n",(double)sum/cnt);
		}

		void buildLockGraphRecur(int n1,int n2,vector<RecurseInfo> &info,int depth=0) {
			if (depth >= maxDepth) {
				info.resize(info.size()+1);
				info[info.size()-1].n1=n1;
				info[info.size()-1].n2=n2;
				return;
			}
			if(n1==n2) {
				KDNode &p1=pool[n1];
				bool lock1=p1.numParts>=0 || p1.size<2.0*(p1.maxRad+p1.searchRadius);
				if(!lock1) {
					buildLockGraphRecur(p1.firstChild,p1.firstChild,info,depth+1);
					buildLockGraphRecur(p1.firstChild,p1.firstChild+1,info,depth+1);
					buildLockGraphRecur(p1.firstChild+1,p1.firstChild+1,info,depth+1);
				}
			} else {
				KDNode &p1=pool[n1];
				KDNode &p2=pool[n2];
				double dx=0.5*(p1.min[0]+p1.max[0]-p2.min[0]-p2.max[0]);
				double dy=0.5*(p1.min[1]+p1.max[1]-p2.min[1]-p2.max[1]);
				double dz=0.5*(p1.min[2]+p1.max[2]-p2.min[2]-p2.max[2]);
				double dsqr=dx*dx+dy*dy+dz*dz;
				double dist=sqrt(dsqr);
				double sdist=p1.searchRadius+p1.maxRad+p2.searchRadius+p2.maxRad;
				if(dist-0.866*(p1.size+p2.size)<=sdist) {
					bool lock1=p1.numParts>=0 || p1.size<2.0*(p1.maxRad+p1.searchRadius);
					bool lock2=p2.numParts>=0 || p2.size<2.0*(p2.maxRad+p2.searchRadius);
					if(lock1 && lock2) {
						addLockEdge(n1,n2);
					} else if(lock1) {
						buildLockGraphRecur(n1,p2.firstChild,info,depth+1);
						buildLockGraphRecur(n1,p2.firstChild+1,info,depth+1);
					} else if(lock2) {
						buildLockGraphRecur(n2,p1.firstChild,info,depth+1);
						buildLockGraphRecur(n2,p1.firstChild+1,info,depth+1);
					} else {
						buildLockGraphRecur(p1.firstChild,p2.firstChild,info,depth+1);
						buildLockGraphRecur(p1.firstChild,p2.firstChild+1,info,depth+1);
						buildLockGraphRecur(p1.firstChild+1,p2.firstChild,info,depth+1);
						buildLockGraphRecur(p1.firstChild+1,p2.firstChild+1,info,depth+1);
					}
				}
			}
		}

		void addLockEdge(int n1,int n2) {
			omp_set_lock(&lock);
			if(nodeEdges[n1].size()==0) {
				setParticleLocks(n1,n1);
			}
			if(nodeEdges[n2].size()==0) {
				setParticleLocks(n2,n2);
			}
			nodeEdges[n1].push_back(n2);
			nodeEdges[n2].push_back(n1);
			omp_unset_lock(&lock);
		}

		void setParticleLocks(int n,int lockNode) {
			if(pool[n].numParts>=0) {
				for(int i=0; i<pool[n].numParts; ++i) {
					particleLockNodes[pool[n].parts[i].i]=lockNode;
				}
			} else {
				setParticleLocks(pool[n].firstChild,lockNode);
				setParticleLocks(pool[n].firstChild+1,lockNode);
			}
		}
#endif

		template<class Population>
		void buildRecur(Population &pop,vector<ParticleIndex> &index,int start,int end,int node,vector<RecurseInfo> &info,int depth=0) {
			if (depth >= maxDepth) {
				info.resize(info.size()+1);
				info[info.size()-1].n1=node;
				info[info.size()-1].n2=start;
				info[info.size()-1].end=end;
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
				double min[]={pop.getx(index[start]),pop.gety(index[start]),pop.getz(index[start])};
				double max[]={pop.getx(index[start]),pop.gety(index[start]),pop.getz(index[start])};
				for(int i=start+1; i<end; ++i) {
					for(int j=0; j<DIM; ++j) {
						double v=pop.get(index[i],j);
						if(v<min[j]) min[j]=v;
						if(v>max[j]) max[j]=v;
					}
				}
				int dim=0;
				for(int i=1; i<DIM; ++i) {
					if(max[i]-min[i]>max[dim]-min[dim]) dim=i;
				}
				//printf("Split %d %d in %d\n",end,start,dim);
				int mid=(start+end)/2;
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
				//printf("Node %d %d : %e = %e-%e, %e-%e, %e-%e\n",firstFree,dim,pop.get(index[mid],dim),max[0],min[0],max[1],min[1],max[2],min[2]);
				buildRecur(pop,index,start,mid+1,pool[node].firstChild,info,depth+1);
				buildRecur(pop,index,mid+1,end,pool[node].firstChild+1,info,depth+1);
			}
		}

		template<class Population>
		void findMid(Population &pop,vector<ParticleIndex> &index,int start,int end,int mid,int dim) {
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
					if(pool[fc].min[i]<pool[n].min[i])
						printf("Invalid min0 node=%d child=%d dim=%n",n,fc,i);
					if(pool[fc].max[i]>pool[n].max[i])
						printf("Invalid min0 node=%d child=%d dim=%n",n,fc,i);
					if(pool[fc+1].min[i]<pool[n].min[i])
						printf("Invalid min1 node=%d child=%d dim=%n",n,fc+1,i);
					if(pool[fc+1].max[i]>pool[n].max[i])
						printf("Invalid min1 node=%d child=%d dim=%n",n,fc+1,i);
				}
				if(pool[fc].min[pool[n].splitDim]>pool[n].splitVal)
					printf("Invalid side of split min0 node=%d child=%d dim=%d\n",n,fc,pool[n].splitDim);
				if(pool[fc].max[pool[n].splitDim]>pool[n].splitVal)
					printf("Invalid side of split max0 node=%d child=%d dim=%d\n",n,fc,pool[n].splitDim);
				if(pool[fc+1].min[pool[n].splitDim]<pool[n].splitVal)
					printf("Invalid side of split min1 node=%d child=%d dim=%d\n",n,fc+1,pool[n].splitDim);
				if(pool[fc+1].max[pool[n].splitDim]<pool[n].splitVal)
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
					if(pop.getx(pool[n].parts[i])<pool[n].min[0])
						printf("Invalid x node=%d part=%d\n",n,pool[n].parts[i]);
					if(pop.getx(pool[n].parts[i])>pool[n].max[0])
						printf("Invalid x node=%d part=%d\n",n,pool[n].parts[i]);
					if(pop.gety(pool[n].parts[i])<pool[n].min[1])
						printf("Invalid y node=%d part=%d\n",n,pool[n].parts[i]);
					if(pop.gety(pool[n].parts[i])>pool[n].max[1])
						printf("Invalid y node=%d part=%d\n",n,pool[n].parts[i]);
					if(pop.getz(pool[n].parts[i])<pool[n].min[2])
						printf("Invalid z node=%d part=%d\n",n,pool[n].parts[i]);
					if(pop.getz(pool[n].parts[i])>pool[n].max[2])
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
				for(int i=0; i<DIM; ++i) {
					pool[n].min[i]=min(pool[pool[n].firstChild].min[i],pool[pool[n].firstChild+1].min[i]);
					pool[n].max[i]=max(pool[pool[n].firstChild].max[i],pool[pool[n].firstChild+1].max[i]);
				}
				pool[n].maxRad=max(pool[pool[n].firstChild].maxRad,pool[pool[n].firstChild+1].maxRad);
			} else if(pool[n].numParts>0) {
				pool[n].maxRad=0.0;
				int num=pool[n].numParts;
				for(int i=0; i<DIM; ++i) {
					pool[n].min[i]=pop.get(pool[n].parts[0],i);
					pool[n].max[i]=pop.get(pool[n].parts[0],i);
				}
				pool[n].velDispSum=0.0;
				pool[n].cnt=0.0;
				for(int i=0; i<num; ++i) {
					double radius=pop.getRadius(pool[n].parts[i]);
					if(radius>pool[n].maxRad) pool[n].maxRad=radius;
					for(int j=i+1; j<num; ++j) {
						double dx=pop.getvx(pool[n].parts[i])-pop.getvx(pool[n].parts[j]);
						double dy=pop.getvy(pool[n].parts[i])-pop.getvy(pool[n].parts[j]);
						double dz=pop.getvz(pool[n].parts[i])-pop.getvz(pool[n].parts[j]);
						pool[n].velDispSum+=dx*dx+dy*dy+dz*dz;
						pool[n].cnt+=1.0;
//						printf("%e %e %e %e %e\n",dx,dy,dz,searchRadius,cnt);
					}
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
						if(x[k]<pool[n].min[k]) pool[n].min[k]=x[k];
						if(x[k]>pool[n].max[k]) pool[n].max[k]=x[k];
					}
				}
			} else {
				pool[n].maxRad=0.0;
				pool[n].velDispSum=0.0;
				pool[n].cnt=0.0;
			}
			if(pool[n].cnt>0) pool[n].searchRadius=5*sqrt(pool[n].velDispSum/pool[n].cnt)*pop.getTimeStep();
			else pool[n].searchRadius=0.0;
			double diff=pool[n].max[0]-pool[n].min[0];
			for(int i=1; i<DIM; ++i) {
				double d=pool[n].max[i]-pool[n].min[i];
				if(d>diff) diff=d;
			}
			pool[n].size=diff;
//			printf("Node %d = %e-%e, %e-%e, %e-%e, %e\n",n,pool[n].min[0],pool[n].max[0],pool[n].min[1],pool[n].max[1],pool[n].min[2],pool[n].max[2],pool[n].size);
#ifdef GRAVITY
			if(pool[n].mass>0.0) {
				pool[n].cmx/=pool[n].mass;
				pool[n].cmy/=pool[n].mass;
				pool[n].cmz/=pool[n].mass;
			}
#endif
		}

		template<class Population,typename ActionType>
		void recurseForAll(Population &pop,ActionType &action,int node1,int node2,std::vector<RecurseInfo> &info,int depth=0) {
			if (depth >= maxDepth)
			{
				info.resize(info.size()+1);
				info[info.size()-1].n1=node1;
				info[info.size()-1].n2=node2;
				return;
			}
			//printf("%d %d %e\n b1=(%e %e, %e %e)\n b2=(%e %e, %e %e)\n",node1,node2,searchRadius,b1.min[0],b1.max[0],b1.min[1],b1.max[1],b2.min[0],b2.max[0],b2.min[1],b2.max[1]);
			if(pool[node1].numParts>=0 && pool[node2].numParts>=0) {
				int num1=pool[node1].numParts;
				for(int i=0; i<num1; ++i) {
					int j;
					if(node1==node2) j=i+1;
						else j=0;
					int num2=pool[node2].numParts;
					for(; j<num2; ++j) {
						action(pop,pool[node1].parts[i],pool[node2].parts[j]);
					}
				}
			} else if(pool[node1].numParts>=0) {
				recurseForAll(pop,action,node2,node1,info,depth+1);
			} else if(node1==node2) {
				// This is a special case that prevents me from
				// doing extra work. In this situation I know they
				// all overlap and I only need to recurse in 3 different
				// combination, not 4.
				recurseForAll(pop,action,pool[node1].firstChild,pool[node2].firstChild,info,depth+1);
				recurseForAll(pop,action,pool[node1].firstChild+1,pool[node2].firstChild+1,info,depth+1);
				recurseForAll(pop,action,pool[node1].firstChild,pool[node2].firstChild+1,info,depth+1);
			} else {
				int sd=pool[node1].splitDim;
				double sv=pool[node1].splitVal;
				if(pool[node2].min[sd]-(max(pool[node1].searchRadius,pool[node2].searchRadius)+pool[node1].maxRad+pool[node2].maxRad)<sv) {
					recurseForAll(pop,action,node2,pool[node1].firstChild,info,depth+1);
				}
				if(pool[node2].max[sd]+(max(pool[node1].searchRadius,pool[node2].searchRadius)+pool[node1].maxRad+pool[node2].maxRad)>sv) {
					recurseForAll(pop,action,node2,pool[node1].firstChild+1,info,depth+1);
				}
			}
		}

		template<class Population,class CollisionManager>
		void recurseForOne(Population &pop,CollisionManager &cm,ParticleIndex p,ParticleIndex ignore,int node) {
			if(pool[node].numParts>=0) {
				int num=pool[node].numParts;
				for(int i=0; i<num; ++i) {
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
				double val = location[p.i].c[sd];
				if(val-(pool[node].searchRadius+pool[node].maxRad+pop.getRadius(p))<sv) {
					recurseForOne(pop,cm,p,ignore,pool[node].firstChild);
				}
				if(val+(pool[node].searchRadius+pool[node].maxRad+pop.getRadius(p))>sv){
					recurseForOne(pop,cm,p,ignore,pool[node].firstChild+1);
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
		void doForce(int n, Population &pop, ParticleIndex pi, double offsetX, double offsetY, vector<AccelVect> &acc) {
			if(pool[n].mass<=0.0) return;
			int num=pool[n].numParts;
			if(num>0) {
				for(int j=0; j<num; ++j) {
					ParticleIndex oi=pool[n].parts[j];
					if(oi!=pi && oi>=0) {
						double dx=pop.getx(oi)+offsetX-pop.getx(pi);
						double dy=pop.gety(oi)+offsetY-pop.gety(pi);
						double dz=pop.getz(oi)-pop.getz(pi);
						double dist=sqrt(dx*dx+dy*dy+dz*dz);
#ifdef GRAV_STATS
						IterCounter::recordSqr();
						IterCounter::recordSqrt();
#endif
						acc[pi.i] += shortRangeForce.applyForce(pop, pi, oi, dx, dy, dz, dist, offsetX, offsetY);
					}
#ifdef PARALLEL
					else if(oi<0) {
//						printf("rp force %d of %d, %e %e %e %e\n",j,pool[n].numParts,remoteParticles[-oi].x,remoteParticles[-oi].y,remoteParticles[-oi].z,remoteParticles[-oi].mass);
						double dx=remoteParticles[-oi.i].x+offsetX-pop.getx(pi);
						double dy=remoteParticles[-oi.i].y+offsetY-pop.gety(pi);
						double dz=remoteParticles[-oi.i].z-pop.getz(pi);
						double dist=sqrt(dx*dx+dy*dy+dz*dz);
#ifdef GRAV_STATS
						IterCounter::recordSqr();
						IterCounter::recordSqrt();
#endif
						if(dist>2*pop.getRadius(pi)) {
							double mag=remoteParticles[-oi.i].mass/(dist*dist*dist);
							acc[pi.i].ax+=dx*mag;
							acc[pi.i].ay+=dy*mag;
							acc[pi.i].az+=dz*mag;
						}
					}
#endif
				}
			} else if(num<0){
				double dx=pool[n].cmx+offsetX-pop.getx(pi);
				double dy=pool[n].cmy+offsetY-pop.gety(pi);
				double dz=pool[n].cmz-pop.getz(pi);
				double dsqr=dx*dx+dy*dy+dz*dz;
#ifdef GRAV_STATS
				IterCounter::recordSqr();
#endif
				if(theta*theta*dsqr>pool[n].size*pool[n].size || pool[n].firstChild<0) {
					double dist=sqrt(dsqr);
#ifdef GRAV_STATS
					IterCounter::recordSqrt();
#endif
					double mag = pool[n].mass/(dist*dsqr);
					acc[pi.i].ax += dx*mag;
					acc[pi.i].ay += dy*mag;
					acc[pi.i].az += dz*mag;
				} else {
					if(dsqr+pool[n].size*pool[n].size>minDist*minDist) {
						doForce(pool[n].firstChild,pop,pi,offsetX,offsetY, acc);
						doForce(pool[n].firstChild+1,pop,pi,offsetX,offsetY, acc);
					}
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
			if(pool[pullerNode].numParts>0) {
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
#ifdef GRAV_STATS
			IterCounter::recordSqr();
			IterCounter::recordSqrt();
#endif
					if(theta*(ndist-pool[pulledNode].size)>pool[pullerNode].size) {
						// Far away so pull all
						if(depth>=0) {
							GravRecurseInfo gri(pullerNode,pulledNode,offsetX,offsetY);
							stack.push_back(gri);
							return;
						}
						pullAllChildren(pullerNode,pulledNode,pop,offsetX,offsetY,acc);
					} else if(theta*(ndist+pool[pulledNode].size)<pool[pullerNode].size || pool[pulledNode].numParts>0) {
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
				int numk=pool[pullerNode].numParts;
				if(numk>0) {
					for(int k=0; k<numk; ++k) {
						ParticleIndex oi=pool[pullerNode].parts[k];
						if(oi>=0 && oi<pop.getNumReal()) {
							int j=0;
							if(pullerNode==pulledNode && offsetX==0.0 && offsetY==0) j=k+1;
							int num=pool[pulledNode].numParts;
							for(; j<num; ++j) {
								ParticleIndex pi=pool[pulledNode].parts[j];
								double dx=pop.getx(oi)+offsetX-pop.getx(pi);
								double dy=pop.gety(oi)+offsetY-pop.gety(pi);
								double dz=pop.getz(oi)-pop.getz(pi);
								double dist=sqrt(dx*dx+dy*dy+dz*dz);
	#ifdef GRAV_STATS
				IterCounter::recordSqr();
				IterCounter::recordSqrt();
	#endif
								double mag=pop.getMass(oi)/(dist*dist*dist);
								acc[pi.i].ax+=dx*mag;
								acc[pi.i].ay+=dy*mag;
								acc[pi.i].az+=dz*mag;
								if(pullerNode==pulledNode) {
									mag=-pop.getMass(pi)/(dist*dist*dist);
									acc[oi.i].ax+=dx*mag;
									acc[oi.i].ay+=dy*mag;
									acc[oi.i].az+=dz*mag;
								}
							}
						}
#ifdef PARALLEL
						else if(oi<0) {
							for(int j=0; j<pool[pulledNode].numParts; ++j) {
								ParticleIndex pi(pool[pulledNode].parts[j]);
								double dx=remoteParticles[-oi.i].x+offsetX-pop.getx(pi);
								double dy=remoteParticles[-oi.i].y+offsetY-pop.gety(pi);
								double dz=remoteParticles[-oi.i].z-pop.getz(pi);
								double dist=sqrt(dx*dx+dy*dy+dz*dz);
								if(dist>2*pop.getRadius(pi)) {
									double mag=remoteParticles[-oi.i].mass/(dist*dist*dist);
									acc[pi.i].ax+=dx*mag;
									acc[pi.i].ay+=dy*mag;
									acc[pi.i].az+=dz*mag;
								}
							}
						}
#endif
					}
				} else {
					int num=pool[pulledNode].numParts;
					for(int k=0; k<num; ++k) {
						ParticleIndex pi=pool[pulledNode].parts[k];
						double dx=pool[pullerNode].cmx+offsetX-pop.getx(pi);
						double dy=pool[pullerNode].cmy+offsetY-pop.gety(pi);
						double dz=pool[pullerNode].cmz-pop.getz(pi);
						double dsqr=dx*dx+dy*dy+dz*dz;
						double dist=sqrt(dsqr);
#ifdef GRAV_STATS
			IterCounter::recordSqr();
			IterCounter::recordSqrt();
#endif
						double mag=pool[pullerNode].mass/(dist*dsqr);
						acc[pi.i].ax+=dx*mag;
						acc[pi.i].ay+=dy*mag;
						acc[pi.i].az+=dz*mag;
					}
				}
			} else {
				pullAllChildren(pullerNode,pool[pulledNode].firstChild,pop,offsetX,offsetY,acc);
				pullAllChildren(pullerNode,pool[pulledNode].firstChild+1,pop,offsetX,offsetY,acc);
			}
		}

		template<class Population>
		void reorderFromTree(Population &pop,int n,vector<int> &map,vector<int> &invmap,int &nextPos) {
			if(pool[n].numParts<0) {
				reorderFromTree(pop,pool[n].firstChild,map,invmap,nextPos);
				reorderFromTree(pop,pool[n].firstChild+1,map,invmap,nextPos);
			} else {
				int num=pool[n].numParts;
				for(int i=0; i<num; ++i) {
					int orig=pool[n].parts[i];
					if(nextPos!=map[orig]) {
						pop.swapParticles(nextPos,map[orig]);
						std::swap(location[nextPos],location[map[orig]]);
						int tmp=invmap[nextPos];
						std::swap(invmap[nextPos],invmap[map[orig]]);
						std::swap(map[tmp],map[orig]);
					}
					pool[n].parts[i]=nextPos;
					nextPos++;
				}
			}
		}

#ifdef PARALLEL
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
						boundsBuffer[2*j]=pool[0].min[j];
						boundsBuffer[2*j+1]=pool[0].max[j];
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
				double xSep=std::max(xmin-(pool[n].max[0]+offsetX),(pool[n].min[0]+offsetX)-xmax);
				if(xSep<0.0) xSep=0.0;
				double ySep=std::max(ymin-(pool[n].max[1]+offsetY),(pool[n].min[1]+offsetY)-ymax);
				if(ySep<0.0) ySep=0.0;
				double zSep=std::max(zmin-(pool[n].max[2]),(pool[n].min[2])-zmax);
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
				treeBuffer.push_back(pool[n].min[0]+offsetX);
				treeBuffer.push_back(pool[n].max[0]+offsetX);
				treeBuffer.push_back(pool[n].min[1]+offsetY);
				treeBuffer.push_back(pool[n].max[1]+offsetY);
				treeBuffer.push_back(pool[n].min[2]);
				treeBuffer.push_back(pool[n].max[2]);
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
					printf("%d - Resize pool from %d\n",pc.getProcessNum(),pool.size());
					pool.resize(pool.size()+pool.size()/2);
					printf("%d - Resized pool to %d\n",pc.getProcessNum(),pool.size());
				}
				pool[firstFree].firstChild=-1;
				if(treeBuffer[i]==0) {
					//printf("CM node\n");
					pool[firstFree].cmx=treeBuffer[i+1];
					pool[firstFree].cmy=treeBuffer[i+2];
					pool[firstFree].cmz=treeBuffer[i+3];
					pool[firstFree].mass=treeBuffer[i+4];
					pool[firstFree].min[0]=treeBuffer[i+5];
					pool[firstFree].max[0]=treeBuffer[i+6];
					pool[firstFree].min[1]=treeBuffer[i+7];
					pool[firstFree].max[1]=treeBuffer[i+8];
					pool[firstFree].min[2]=treeBuffer[i+9];
					pool[firstFree].max[2]=treeBuffer[i+10];
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
							pool[firstFree].min[0]=x;
							pool[firstFree].max[0]=x;
							pool[firstFree].min[1]=y;
							pool[firstFree].max[1]=y;
							pool[firstFree].min[2]=z;
							pool[firstFree].max[2]=z;
						} else {
							if(x<pool[firstFree].min[0]) pool[firstFree].min[0]=x;
							if(x>pool[firstFree].max[0]) pool[firstFree].max[0]=x;
							if(y<pool[firstFree].min[1]) pool[firstFree].min[1]=y;
							if(y>pool[firstFree].max[1]) pool[firstFree].max[1]=y;
							if(z<pool[firstFree].min[2]) pool[firstFree].min[2]=z;
							if(z>pool[firstFree].max[2]) pool[firstFree].max[2]=z;
						}
					}
				} else {
					printf("Error: Bad code in treeBuffer! %d is %f\n",i,treeBuffer[i]);
					exit(-1);
				}
//				printf("Bounds %e %e %e %e\n",
//					pool[firstFree].minx,pool[firstFree].maxx,
//					pool[firstFree].miny,pool[firstFree].maxy);
				pool[firstFree].size=std::max(std::max(pool[firstFree].max[0]-pool[firstFree].min[0],pool[firstFree].max[1]-pool[firstFree].min[1]),pool[firstFree].max[2]-pool[firstFree].min[2]);
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
							printf("new bounds %e %e %e %e %e %e\n",pool[firstFree].min[0],pool[firstFree].max[0],pool[firstFree].min[1],pool[firstFree].max[1],pool[firstFree].min[2],pool[firstFree].max[2]);
							printf("root bounds %e %e %e %e %e %e\n",pool[curRoot].min[0],pool[curRoot].max[0],pool[curRoot].min[1],pool[curRoot].max[1],pool[curRoot].min[2],pool[curRoot].max[2]);
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
			return pool[n].min[0]>=pool[m].min[0] &&
				pool[n].max[0]<=pool[m].max[0] &&
				pool[n].min[1]>=pool[m].min[1] &&
				pool[n].max[1]<=pool[m].max[1] &&
				pool[n].min[2]>=pool[m].min[2] &&
				pool[n].max[2]<=pool[m].max[2];
		}

		void addNode(int n,int newNode) {
			int p=-1;
			double cx=(pool[newNode].max[0]+pool[newNode].min[0])*0.5;
			double cy=(pool[newNode].max[1]+pool[newNode].min[1])*0.5;
			double cz=(pool[newNode].max[2]+pool[newNode].min[2])*0.5;
			while(n!=-1) {
				//printf("%d %d %d\n",pc.getProcessNum(),n,p);
				p=n;
				n=pool[n].childNode(cx,cy,cz);
				if(n==0) {
					printf("Error %d: add node got to 0, p=%d, newNode=%d, pool.size=%d.\n",pc.getProcessNum(),p,newNode,pool.size());
					printf("p.numParts=%d, n.numParts=%d\n",pool[p].numParts,pool[newNode].numParts);
					printf("cx=%e, cy=%e\n",cx,cy);
					printf("p.splitDim=%d, p.splitVal=%e\n",pool[p].splitDim,pool[p].splitVal);
					printf("p bounds %e %e %e %e %e %e\n",pool[p].min[0],pool[p].max[0],pool[p].min[1],pool[p].max[1],pool[p].min[2],pool[p].max[2]);
					printf("new bounds %e %e %e %e %e %e\n",pool[newNode].min[0],pool[newNode].max[0],pool[newNode].min[1],pool[newNode].max[1],pool[newNode].min[2],pool[newNode].max[2]);
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
				printf("p bounds %e %e %e %e\n",pool[p].min[0],pool[p].max[0],pool[p].min[1],pool[p].max[1]);
				printf("new bounds %e %e %e %e\n",pool[newNode].min[0],pool[newNode].max[0],pool[newNode].min[1],pool[newNode].max[1]);
			}
		}
#endif

		vector<KDNode> pool;
		vector<Location> location;
#ifndef LOCK_GRAPH
		vector<vector<bool> > inUse;
#else
		vector<vector<int> > nodeEdges;
		vector<int> particleLockNodes;
#endif
		vector<vector<AccelVect> > acc;
		double gridSpacing;
		int maxDepth;
		int firstFree;
		double theta;
		Boundary &bounds;  // boundary conditions
		ShortRangeForce shortRangeForce;
		double minDist;
		omp_lock_t lock;
#ifdef REORDER
		int stepCount;
#endif
#ifdef PARALLEL
		ProcessorCommunication &pc;
		vector<double> boundsBuffer;
		vector<double> sizeBuffer;
		vector<double> treeBuffer;
		vector<int> receivedTrees;
		vector<RemoteParticle> remoteParticles;
#endif
};

#endif