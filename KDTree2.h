// This is my version of a spatial KD tree used for collision detection.

#include <vector>
#include <algorithm>

#define MAX_NUM 2
#define DIM 3

using std::min;
using std::max;

struct KDTreeNode {
	int particle[MAX_NUM];
	int numParts;
	int splitDim;
	double splitValue;
	int left,right;
	double maxRad;
};

#ifndef PARALLEL_GRAV_COLL_TREE
struct KDBounds {
	double min[DIM],max[DIM];
};

struct Location {
	double c[DIM];
};

struct RecurseInfo {
	int n1,n2;
	KDBounds b1,b2;
};
#endif

class KDTree {
	public:
		KDTree() {
			maxDepth = (int) (log(omp_get_max_threads())/log(2) + 2);
		}
		template<class Population>
		void build(Population &pop) {
//			printf("Start build\n");
			pool.resize(pop.getNumBodies()*2);
			location.resize(pop.getNumBodies());
			root=0;
			pool[root].particle[0]=0;
			pool[root].numParts=1;
			for (int i = 0; i < DIM; i++)
				location[0].c[i] = pop.get(0,i);
			firstFree=1;
			for(int i=0; i<DIM; ++i) {
				totalBounds.min[i]=1e100;
				totalBounds.max[i]=-1e100;
			}

			for(int i=1; i<pop.getNumBodies(); ++i) {
				add(pop,i);
				for(int j=0; j<DIM; ++j) {
					location[i].c[j] = pop.get(i,j);
					if (pop.get(i,j) < totalBounds.min[j])
						totalBounds.min[j]=pop.get(i,j);
					if (pop.get(i,j) > totalBounds.max[j])
						totalBounds.max[j]=pop.get(i,j);
				}
			}

			finalize(pop,root);
			printf("End build %d nodes\n",firstFree);
			fflush(stderr);
//			printTree(pop,root);
		}

		double getMinSpacing() { return 3*searchRadius; }

		template<class Population,class CollisionManager>
		void findAllCollisions(Population &pop,CollisionManager &cm) {
			KDBounds b1,b2;
			for(int i=0; i<DIM; ++i) {
				b1.min[i]=totalBounds.min[i];
				b1.max[i]=totalBounds.max[i];
				b2.min[i]=totalBounds.min[i];
				b2.max[i]=totalBounds.max[i];
			}		
//			printf("Start all\n");
			cnt=0.0;
			searchRadius=0.0;
			setSearchRadius(pop,root);
			printf("%e %e\n",cnt,searchRadius);
			searchRadius=5.0*sqrt(searchRadius/cnt)*pop.getTimeStep();
			printf("%e\n",searchRadius);
			buildLockGrid(pop);

			//Recurse through the top maxDepth levels of the tree.
			//Instead of making further recursive calls, store the parameters for
			//each call as a RecurseInfo and store it in the info queue.

			std::vector<RecurseInfo> info(0);
			recurseForTop(pop,cm,root,b1,root,b2,info);
			
			//Make all stored recursive calls in parallel.

			int isize = (int) info.size();
			#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < isize; i++)
				recurseForAll(pop,cm,info[i].n1,info[i].b1,info[i].n2,info[i].b2);
//			printf("End all\n");
		}

		template<class Population,class CollisionManager>
		void findCollisions(Population &pop,CollisionManager &cm,int p,int ignore) {
			KDBounds b;
			for(int i=0; i<DIM; ++i) {
				b.min[i]=totalBounds.min[i];
				b.max[i]=totalBounds.max[i];
			}			
			recurseForOne(pop,cm,p,ignore,b,root);
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
			gridSpacing = 2*pool[root].maxRad + searchRadius;
			int xlen = (int) ((totalBounds.max[0] - totalBounds.min[0])/gridSpacing + 1);
			int ylen = (int) ((totalBounds.max[1] - totalBounds.min[1])/gridSpacing + 1);
			inUse.resize(xlen);
			for (int i = 0; i < xlen; i++)
			{
				inUse[i].resize(ylen);
				for (int j = 0; j < ylen; j++)
					inUse[i][j] = false;
			}
		}

		template<class Population>
		void add(Population &pop,int p) {
			int node=root;
//			printf("Find node to put in\n");
			while(pool[node].numParts==-1) {
				double val=pop.get(p,pool[node].splitDim);
				if(val<pool[node].splitValue) {
					node=pool[node].left;
				} else {
					node=pool[node].right;
				}
			}
//			printf("Adding to node %d, numparts is %d\n",node,pool[node].numParts);
			if(pool[node].numParts>=MAX_NUM) {
				double min[DIM];
				double max[DIM];
				for(int i=0; i<DIM; ++i) {
					min[i]=pop.get(p,i);
					max[i]=min[i];
				}
				for(int i=0; i<MAX_NUM; ++i) {
					for(int j=0; j<DIM; ++j) {
						double val=pop.get(pool[node].particle[i],j);
						if(val<min[j]) min[j]=val;
						if(val>max[j]) max[j]=val;
					}
				}
				int sd=0;
				double sv;
				for(int i=1; i<DIM; ++i) {
					if(max[i]-min[i]>max[sd]-min[sd]) sd=i;
				}
				sv=0.5*(max[sd]+min[sd]);
				pool[node].splitDim=sd;
				pool[node].splitValue=sv;
				pool[node].left=firstFree;
				pool[node].right=firstFree+1;
				pool[firstFree].numParts=0;
				pool[firstFree+1].numParts=0;
				for(int i=0; i<MAX_NUM; ++i) {
					double val=pop.get(pool[node].particle[i],sd);
					int dest=firstFree;
					if(val>sv) dest=firstFree+1;
					pool[dest].particle[pool[dest].numParts]=
						pool[node].particle[i];
					pool[dest].numParts++;
				}
				pool[node].numParts=-1;
				double val=pop.get(p,sd);
				node=firstFree;
				if(val>sv) node=firstFree+1;
				firstFree+=2;
//				printf("After split %d %d\n",node,firstFree);
			}
			pool[node].particle[pool[node].numParts]=p;
			pool[node].numParts++;
		}

		template<class Population>
		void finalize(Population &pop,int node) {
//			printf("Finalizing %d - %d\n",node,pool[node].numParts);
			if(pool[node].numParts>=0) {
				pool[node].maxRad=0.0;
				for(int i=1; i<pool[node].numParts; ++i) {
					double radius=pop.getRadius(pool[node].particle[i]);
					if(radius>pool[node].maxRad) pool[node].maxRad=radius;
				}				
			} else {
				int left=pool[node].left;
				int right=pool[node].right;
				finalize(pop,left);
				finalize(pop,right);
				pool[node].maxRad=max(pool[left].maxRad,pool[right].maxRad);
			}
		}

		template<class Population,class CollisionManager>
		void recurseForTop(Population &pop, CollisionManager &cm, int node1, KDBounds &b1, int node2, KDBounds &b2, std::vector<RecurseInfo> &info, int depth=0) {
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
						double t=pop.collisionTime(pool[node1].particle[i],pool[node2].particle[j]);
						if(t>=0.0 && t<pop.getTimeStep()) {
							cm.addPotentialCollision(pool[node1].particle[i],pool[node2].particle[j],t);
						}
					}
				}
			} else if(pool[node1].numParts>=0) {
				recurseForTop(pop,cm,node2,b2,node1,b1,info,depth+1);
			} else if(node1==node2) {
				// This is a special case that prevents me from
				// doing extra work. In this situation I know they
				// all overlap and I only need to recurse in 3 different
				// combination, not 4.
				int sd=pool[node1].splitDim;
				double sv=pool[node1].splitValue;
				double maxTmp=b1.max[sd];
				double minTmp=b1.min[sd];
				
				b1.max[sd]=b2.max[sd]=sv;
				recurseForTop(pop,cm,pool[node1].left,b1,pool[node2].left,b2,info,depth+1);
				b1.max[sd]=b2.max[sd]=maxTmp;
				b1.min[sd]=b2.min[sd]=sv;
				recurseForTop(pop,cm,pool[node1].right,b1,pool[node2].right,b2,info,depth+1);
				b1.min[sd]=minTmp;
				b1.max[sd]=sv;
				recurseForTop(pop,cm,pool[node1].left,b1,pool[node2].right,b2,info,depth+1);
				b1.max[sd]=maxTmp;
				b2.min[sd]=minTmp;
			} else {
				int sd=pool[node1].splitDim;
				double sv=pool[node1].splitValue;
				double tmp;
//				int pick=0;
				if(b2.min[sd]-(searchRadius+pool[node1].maxRad+pool[node2].maxRad)<sv) {
					tmp=b1.max[sd];
					b1.max[sd]=sv;
					recurseForTop(pop,cm,node2,b2,pool[node1].left,b1,info,depth+1);
					b1.max[sd]=tmp;
//					pick|=1;
				}
				if(b2.max[sd]+(searchRadius+pool[node1].maxRad+pool[node2].maxRad)>sv) {
					tmp=b1.min[sd];
					b1.min[sd]=sv;
					recurseForTop(pop,cm,node2,b2,pool[node1].right,b1,info,depth+1);
					b1.min[sd]=tmp;
//					pick|=2;
				}
//				printf("%d\n",pick);
			}
		}

		template<class Population,class CollisionManager>
		void recurseForAll(Population &pop,CollisionManager &cm,int node1,KDBounds &b1,int node2,KDBounds &b2) {
			//printf("%d %d %e\n b1=(%e %e, %e %e)\n b2=(%e %e, %e %e)\n",node1,node2,searchRadius,b1.min[0],b1.max[0],b1.min[1],b1.max[1],b2.min[0],b2.max[0],b2.min[1],b2.max[1]);
			if(pool[node1].numParts>=0 && pool[node2].numParts>=0) {
				for(int i=0; i<pool[node1].numParts; ++i) {
					int j;
					if(node1==node2) j=i+1;
						else j=0;
					for(; j<pool[node2].numParts; ++j) {
						double t=pop.collisionTime(pool[node1].particle[i],pool[node2].particle[j]);
						if(t>=0.0 && t<pop.getTimeStep()) {
							cm.addPotentialCollision(pool[node1].particle[i],pool[node2].particle[j],t);
						}
					}
				}
			} else if(pool[node1].numParts>=0) {
				recurseForAll(pop,cm,node2,b2,node1,b1);
			} else if(node1==node2) {
				// This is a special case that prevents me from
				// doing extra work. In this situation I know they
				// all overlap and I only need to recurse in 3 different
				// combination, not 4.
				int sd=pool[node1].splitDim;
				double sv=pool[node1].splitValue;
				double maxTmp=b1.max[sd];
				double minTmp=b1.min[sd];
				
				b1.max[sd]=b2.max[sd]=sv;
				recurseForAll(pop,cm,pool[node1].left,b1,pool[node2].left,b2);
				b1.max[sd]=b2.max[sd]=maxTmp;
				b1.min[sd]=b2.min[sd]=sv;
				recurseForAll(pop,cm,pool[node1].right,b1,pool[node2].right,b2);
				b1.min[sd]=minTmp;
				b1.max[sd]=sv;
				recurseForAll(pop,cm,pool[node1].left,b1,pool[node2].right,b2);
				b1.max[sd]=maxTmp;
				b2.min[sd]=minTmp;
			} else {
				int sd=pool[node1].splitDim;
				double sv=pool[node1].splitValue;
				double tmp;
//				int pick=0;
				if(b2.min[sd]-(searchRadius+pool[node1].maxRad+pool[node2].maxRad)<sv) {
					tmp=b1.max[sd];
					b1.max[sd]=sv;
					recurseForAll(pop,cm,node2,b2,pool[node1].left,b1);
					b1.max[sd]=tmp;
//					pick|=1;
				}
				if(b2.max[sd]+(searchRadius+pool[node1].maxRad+pool[node2].maxRad)>sv) {
					tmp=b1.min[sd];
					b1.min[sd]=sv;
					recurseForAll(pop,cm,node2,b2,pool[node1].right,b1);
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
					if(pool[node].particle[i]!=p && pool[node].particle[i]!=ignore) {
						double t=pop.collisionTime(p,pool[node].particle[i]);
						if(t>=0.0 && t<pop.getTimeStep()) {
							cm.addPotentialCollision(p,pool[node].particle[i],t);
						}
					}
				}
			} else {
				int sd=pool[node].splitDim;
				double sv=pool[node].splitValue;
				//double val=pop.get(p,sd);
				double val = location[p].c[sd];
				if(val-(searchRadius+pool[node].maxRad+pop.getRadius(p))<sv) {
					double tmp=b.max[sd];
					b.max[sd]=sv;
					recurseForOne(pop,cm,p,ignore,b,pool[node].left);
					b.max[sd]=tmp;
				}
				if(val+(searchRadius+pool[node].maxRad+pop.getRadius(p))>sv){
					double tmp=b.min[sd];
					b.min[sd]=sv;
					recurseForOne(pop,cm,p,ignore,b,pool[node].right);
					b.min[sd]=tmp;
				}
			}
		}
		
		template<class Population>
		void setSearchRadius(Population &pop,int node) {
			if(pool[node].numParts>=0) {
//				printf("Node %d has %d particles\n",node,pool[node].numParts);
				for(int i=0; i<pool[node].numParts; ++i) {
					for(int j=i+1; j<pool[node].numParts; ++j) {
						double dx=pop.getvx(i)-pop.getvx(j);
						double dy=pop.getvy(i)-pop.getvy(j);
						double dz=pop.getvz(i)-pop.getvz(j);
						searchRadius+=dx*dx+dy*dy+dz*dz;
						cnt+=1.0;
//						printf("%e %e %e %e %e\n",dx,dy,dz,searchRadius,cnt);
					}
				}
			} else {
				setSearchRadius(pop,pool[node].left);
				setSearchRadius(pop,pool[node].right);
			}
		}
		
		template<class Population>
		void printTree(Population &pop,int node) {
			if(pool[node].numParts<0) {
				printf("Node %d %d %e\n",node,pool[node].splitDim,pool[node].splitValue);
				printf("left\n");
				printTree(pop,pool[node].left);
				printf("right\n");
				printTree(pop,pool[node].right);
				printf("done with %d\n",node);
			} else {
				printf("Particles in %d\n",node);
				for(int i=0; i<pool[node].numParts; ++i) {
					int p=pool[node].particle[i];
					printf("%d %e %e %e %e %e %e\n",p,pop.get(p,0),pop.get(p,1)
						,pop.get(p,2),pop.get(p,3),pop.get(p,4),pop.get(p,5));
				}
				printf("Done with %d\n",node);
			}
		}

		std::vector<KDTreeNode> pool;
		std::vector<Location> location;
		std::vector<std::vector<bool> > inUse;
		double gridSpacing;
		int maxDepth;
		int firstFree;
		int root;
		KDBounds totalBounds;
		double searchRadius;
		double cnt;
};
