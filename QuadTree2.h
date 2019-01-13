// This is my version of a spatial KD tree used for collision detection.

#include<vector>
#include<algorithm>

#define MAX_NUM 2

struct QuadTreeNode {
	int particle[MAX_NUM];
	int numParts;
	double cx,cy;
	int child[4];
	double size;
};

class QuadTree {
	public:
		template<class Population>
		void build(Population &pop) {
//			printf("Start build\n");
			pool.resize(pop.getNumBodies()*2);
			root=0;
			pool[root].particle[0]=0;
			pool[root].numParts=1;
			pool[root].child[0]=-1;
			pool[root].cx=(pop.getMaxx()+pop.getMinx())*0.5;
			pool[root].cy=(pop.getMaxy()+pop.getMiny())*0.5;
			pool[root].size=std::max(pop.getMaxx()-pop.getMinx(),pop.getMaxy()-pop.getMiny());
			firstFree=1;
			for(int i=1; i<pop.getNumBodies(); ++i) {
				add(pop,i);
			}
//			printf("End build\n");
//			printTree(pop,root);
		}

		template<class Population,class CollisionManager>
		void findAllCollisions(Population &pop,CollisionManager &cm) {
//			printf("Start all\n");
			cnt=0.0;
			searchRadius=0.0;
			setSearchRadius(pop,root);
//			printf("%e %e\n",cnt,searchRadius);
			searchRadius=2.0*pop.getMaxParticleRadius()+3.0*sqrt(searchRadius/cnt)*pop.getTimeStep();
//			printf("%e\n",searchRadius);
			recurseForAll(pop,cm,root,root);
//			printf("End all\n");
		}

		template<class Population,class CollisionManager>
		void findCollisions(Population &pop,CollisionManager &cm,int p,int ignore) {
//			printf("Find for %d ignore %d\n",p,ignore);
			recurseForOne(pop,cm,p,ignore,root);
		}
	private:
		int child(double x,double y,double cx,double cy) {
			return ((x<cx)?0:1)|((y<cy)?0:2);
		}
		
		template<class Population>
		void add(Population &pop,int p) {
//			printf("Find node to put in\n");
			double x=pop.getx(p);
			double y=pop.gety(p);
			addRecur(pop,p,x,y,root);
		}
		
		template<class Population>
		void addRecur(Population &pop,int p,double x,double y,int node) {
//			printf("Add %d in %d at %e %e\n",p,node,x,y);
			if(pool[node].child[0]>-1) {
				addRecur(pop,p,x,y,pool[node].child[child(x,y,pool[node].cx,pool[node].cy)]);
				return;
			}
//			printf("Adding to node %d, numparts is %d\n",node,pool[node].numParts);
			if(pool[node].numParts<MAX_NUM) {
				pool[node].particle[pool[node].numParts]=p;
				pool[node].numParts++;
			} else {
				pool[node].child[0]=firstFree;
				pool[node].child[1]=firstFree+1;
				pool[node].child[2]=firstFree+2;
				pool[node].child[3]=firstFree+3;
				pool[firstFree].numParts=0;
				pool[firstFree+1].numParts=0;
				pool[firstFree+2].numParts=0;
				pool[firstFree+3].numParts=0;
				pool[firstFree].child[0]=-1;
				pool[firstFree+1].child[0]=-1;
				pool[firstFree+2].child[0]=-1;
				pool[firstFree+3].child[0]=-1;
				double s=pool[node].size*0.5;
				double hs=0.5*s;
				pool[firstFree].cx=pool[node].cx-hs;
				pool[firstFree].cy=pool[node].cy-hs;
				pool[firstFree].size=s;
				pool[firstFree+1].cx=pool[node].cx+hs;
				pool[firstFree+1].cy=pool[node].cy-hs;
				pool[firstFree+1].size=s;
				pool[firstFree+2].cx=pool[node].cx-hs;
				pool[firstFree+2].cy=pool[node].cy+hs;
				pool[firstFree+2].size=s;
				pool[firstFree+3].cx=pool[node].cx+hs;
				pool[firstFree+3].cy=pool[node].cy+hs;
				pool[firstFree+3].size=s;
				firstFree+=4;
//				printf("After split %d %d\n",node,firstFree);
				for(int i=0; i<MAX_NUM; ++i) {
					int pi=pool[node].particle[i];
					double px=pop.getx(pi);
					double py=pop.gety(pi);
					addRecur(pop,pi,px,py,pool[node].child[child(px,py,pool[node].cx,pool[node].cy)]);
				}
				pool[node].numParts=0;
				addRecur(pop,p,x,y,pool[node].child[child(x,y,pool[node].cx,pool[node].cy)]);
			}
		}

		template<class Population,class CollisionManager>
		void recurseForAll(Population &pop,CollisionManager &cm,int node1,int node2) {
			if(pool[node1].child[0]==-1 && pool[node2].child[0]==-1) {
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
			} else if(pool[node1].child[0]==-1) {
				recurseForAll(pop,cm,node2,node1);
			} else if(node1==node2) {
				// This is a special case that prevents me from
				// doing extra work. In this situation I know they
				// all overlap and I only need to recurse in 10 different
				// combination, not 16.
				for(int i=0; i<4; ++i) {
					for(int j=i; j<4; ++j) {
						recurseForAll(pop,cm,pool[node1].child[i],pool[node2].child[j]);
					}
				}
			} else {
				double x1=pool[node1].cx;
				double y1=pool[node1].cy;
				double x2=pool[node2].cx;
				double y2=pool[node2].cy;
				double s2=pool[node2].size*0.5;
				if(x2-s2-searchRadius<x1) {
					if(y2-s2-searchRadius<y1) {
						recurseForAll(pop,cm,node2,pool[node1].child[0]);
					}
					if(y2+s2+searchRadius>y1) {
						recurseForAll(pop,cm,node2,pool[node1].child[2]);
					}
				}
				if(x2+s2+searchRadius>x1) {
					if(y2-s2-searchRadius<y1) {
						recurseForAll(pop,cm,node2,pool[node1].child[1]);
					}
					if(y2+s2+searchRadius>y1) {
						recurseForAll(pop,cm,node2,pool[node1].child[3]);
					}
				}
			}
		}

		template<class Population,class CollisionManager>
		void recurseForOne(Population &pop,CollisionManager &cm,int p,int ignore,int node) {
			if(pool[node].child[0]==-1) {
				for(int i=0; i<pool[node].numParts; ++i) {
					if(pool[node].particle[i]!=p && pool[node].particle[i]!=ignore) {
						double t=pop.collisionTime(p,pool[node].particle[i]);
						if(t>=0.0 && t<pop.getTimeStep()) {
							cm.addPotentialCollision(p,pool[node].particle[i],t);
						}
					}
				}
			} else {
				double x1=pool[node].cx;
				double y1=pool[node].cy;
				double px=pop.getx(p);
				double py=pop.gety(p);
				if(px-searchRadius<x1) {
					if(py-searchRadius<y1) {
						recurseForOne(pop,cm,p,ignore,pool[node].child[0]);
					}
					if(py+searchRadius>y1) {
						recurseForOne(pop,cm,p,ignore,pool[node].child[2]);
					}
				}
				if(px+searchRadius>x1){
					if(py-searchRadius<y1) {
						recurseForOne(pop,cm,p,ignore,pool[node].child[1]);
					}
					if(py+searchRadius>y1) {
						recurseForOne(pop,cm,p,ignore,pool[node].child[3]);
					}
				}
			}
		}
		
		template<class Population>
		void setSearchRadius(Population &pop,int node) {
			if(pool[node].child[0]<0) {
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
				setSearchRadius(pop,pool[node].child[0]);
				setSearchRadius(pop,pool[node].child[1]);
				setSearchRadius(pop,pool[node].child[2]);
				setSearchRadius(pop,pool[node].child[3]);
			}
		}
		/*
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
		*/
		std::vector<QuadTreeNode> pool;
		int firstFree;
		int root;
		double searchRadius;
		double cnt;
};
