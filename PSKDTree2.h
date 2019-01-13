/**
 * PSKDTree.h
 * This file defines a templated phase space KD-tree that can be used to
 * find potential collisions between particles in a population.
 * 
 * This version differs from the first version in that I store the
 * bounds in the nodes instead of passing it down and keeping track
 * of it in the recursion.
**/

#include <vector>
#include <algorithm>

#define MAX_NUM 5
#define DIM 3

struct Bounds {
	double min[2*DIM],max[2*DIM];
};

struct TreeNode {
	int particle[MAX_NUM];
	int numParts;
	int splitDim;
	double splitValue;
	int left,right;
	Bounds b;
	double maxRad;
};

class PSKDTree {
	public:
		template<class Population>
		void build(Population &pop) {
			pool.resize(pop.getNumBodies()*2);
			root=0;
			pool[root].particle[0]=0;
			pool[root].numParts=1;
			firstFree=1;
			for(int i=1; i<pop.getNumBodies(); ++i) {
				add(pop,i);
			}
			finalize(pop,root);
//			printTree(pop,root);
		}

		template<class Population,class CollisionManager>
		void findAllCollisions(Population &pop,CollisionManager &cm) {
			recurseForAll(pop,cm,root,root);
		}

		template<class Population,class CollisionManager>
		void findCollisions(Population &pop,CollisionManager &cm,int p,int ignore) {
			recurseForOne(pop,cm,p,ignore,root);
		}
	private:
		template<class Population>
		void add(Population &pop,int p) {
			int node=root;
//			printf("Find node to put in\n");
			while(pool[node].numParts==0) {
				double val=pop.get(p,pool[node].splitDim);
				if(val<pool[node].splitValue) {
					node=pool[node].left;
				} else {
					node=pool[node].right;
				}
			}
//			printf("Adding to node %d, numparts is %d\n",node,pool[node].numParts);
			if(pool[node].numParts>=MAX_NUM) {
				double min[2*DIM];
				double max[2*DIM];
				for(int i=0; i<2*DIM; ++i) {
					min[i]=pop.get(p,i);
					max[i]=min[i];
				}
				for(int i=0; i<MAX_NUM; ++i) {
					for(int j=0; j<2*DIM; ++j) {
						double val=pop.get(pool[node].particle[i],j);
						if(val<min[j]) min[j]=val;
						if(val>max[j]) max[j]=val;
					}
				}
				int sd=0;
				double sv;
				for(int i=1; i<2*DIM; ++i) {
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
				pool[node].numParts=0;
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
			if(pool[node].numParts>0) {
				for(int d=0; d<2*DIM; ++d) {
					pool[node].b.min[d]=pop.get(pool[node].particle[0],d);
					pool[node].b.max[d]=pop.get(pool[node].particle[0],d);
				}
				pool[node].maxRad=0.0;
				for(int i=1; i<pool[node].numParts; ++i) {
					for(int d=0; d<2*DIM; ++d) {
						double val=pop.get(pool[node].particle[i],d);
						pool[node].b.min[d]<?=val;
						pool[node].b.max[d]>?=val;
					}
					double radius=pop.getRadius(pool[node].particle[i]);
					if(radius>pool[node].maxRad) pool[node].maxRad=radius;
				}				
			} else {
				int left=pool[node].left;
				int right=pool[node].right;
				finalize(pop,left);
				finalize(pop,right);
				for(int d=0; d<2*DIM; ++d) {
					pool[node].b.min[d]=std::min(pool[left].b.min[d],pool[right].b.min[d]);
					pool[node].b.max[d]=std::max(pool[left].b.max[d],pool[right].b.max[d]);
				}
				pool[node].maxRad=std::max(pool[left].maxRad,pool[right].maxRad);
			}
		}

		template<class Population,class CollisionManager>
		void recurseForAll(Population &pop,CollisionManager &cm,int node1,int node2) {
			if(pool[node1].numParts>0 && pool[node2].numParts>0) {
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
			} else if(pool[node1].numParts>0) {
				recurseForAll(pop,cm,node2,node1);
			} else if(node1==node2) {
				// This is a special case that prevents me from
				// doing extra work. In this situation I know they
				// all overlap and I only need to recurse in 3 different
				// combination, not 4.
				recurseForAll(pop,cm,pool[node1].left,pool[node2].left);
				recurseForAll(pop,cm,pool[node1].right,pool[node2].right);
				recurseForAll(pop,cm,pool[node1].left,pool[node2].right);
			} else {
				double dt=pop.getTimeStep();
				double maxRad=std::max(pool[node1].maxRad,pool[node2].maxRad);
				int sd=pool[node1].splitDim;
				int dim=(sd<DIM)?sd:sd-DIM;
				if(checkOverlap(pool[node1].b,pool[node2].b,dt,maxRad,dim)) recurseForAll(pop,cm,node2,pool[node1].left);
				if(checkOverlap(pool[node1].b,pool[node2].b,dt,maxRad,dim)) recurseForAll(pop,cm,node2,pool[node1].right);
			}
		}

		template<class Population,class CollisionManager>
		void recurseForOne(Population &pop,CollisionManager &cm,int p,int ignore,int node) {
			if(pool[node].numParts>0) {
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
				bool below=false;
				bool above=false;
				double val=pop.get(p,sd);
				double dt=pop.getTimeStep();
				double maxRad=pool[node].maxRad;
				Bounds &b=pool[node].b;
				if(sd<DIM) {
					if(val<sv+2.0*maxRad) {
						below=true;
						if(val+pop.get(p,sd+DIM)*dt>sv+dt*b.min[sd+DIM]-2.0*maxRad) above=true;
					}
					if(val>sv-2.0*maxRad){
						above=true;
						if(val+pop.get(p,sd+DIM)*dt<sv+dt*b.max[sd+DIM]+2.0*maxRad) below=true;
					}
				} else {
					if(val<sv+2.0*maxRad/dt) {
						below=true;
						if(pop.get(p,sd-DIM)>=b.min[sd-DIM]-2.0*maxRad) above=true;
					}
					if(val>sv-2.0*maxRad/dt) {
						above=true;
						if(pop.get(p,sd-DIM)<=b.max[sd-DIM]+2.0*maxRad) below=true;
					}
				}
				if(below) {
					recurseForOne(pop,cm,p,ignore,pool[node].left);
				}
				if(above) {
					recurseForOne(pop,cm,p,ignore,pool[node].right);
				}
			}
		}
		
		/**
		 * Determines if particles in two different bounds regions could possibly
		 * overlap.
		 */
		bool checkOverlap(Bounds &b1,Bounds &b2,double dt,double maxRad,int i) {
			return !(((b1.max[i]<b2.min[i]-2.0*maxRad) &&
				(b1.max[i]+dt*b1.max[i+DIM]<b2.min[i]+dt*b2.min[i+DIM]-2.0*maxRad)) ||
				((b1.min[i]>b2.max[i]+2.0*maxRad) &&
				(b1.min[i]+dt*b1.min[i+DIM]>b2.max[i]+dt*b2.max[i+DIM]+2.0*maxRad)));
		}

		template<class Population>
		bool checkOverlap(Bounds &b,Population &pop,int p,double dt,double maxRad) {
			bool ret=true;
			for(int i=0; ret && i<DIM; ++i) {
				if(((b.min[i]>pop.get(p,i)+2.0*maxRad) &&
					(b.min[i]+dt*b.min[i+DIM]>pop.get(p,i)+dt*pop.get(p,i+DIM)+2.0*maxRad)) ||
					((b.max[i]<pop.get(p,i)-2.0*maxRad) &&
					(b.max[i]+dt*b.max[i+DIM]<pop.get(p,i)+dt*pop.get(p,i+DIM)-2.0*maxRad))) {
						ret=false;
				}
			}
			return ret;
		}
		
		template<class Population>
		void printTree(Population &pop,int node) {
			if(pool[node].numParts==0) {
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

		std::vector<TreeNode> pool;
		int firstFree;
		int root;
};
