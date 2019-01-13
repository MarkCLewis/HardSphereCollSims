/**
 * This header file includes the code to do a forcing based on a tree in parallel.
**/

#include <vector>
#include <algorithm>

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

class KDGravTree {
	public:
		KDGravTree(double t,double md=0):theta(t),minDist(md) {}
	
		template<class Population>
		void applyForce(Population &pop) {
			pool.resize(2*pop.getNumBodies());
			pool[0].numParts=0;
			firstFree=1;
			for(int i=0; i<pop.getNumBodies(); ++i) {
				addParticle(0,pop,i);
			}
			finalize(0,pop);
			for(int i=0; i<pop.getNumBodies(); ++i) {
				doForce(0,pop,i);
				pop.adjustAfterForce(i);
			}
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
				printf("Negative node in force. %d\n",n);
				exit(-1);
			}
			if(pool[n].mass<=0.0) return;
//			printf("%d %d %d\n",pc.getProcessNum(),n,i);
			if(pool[n].numParts>0) {
				for(int j=0; j<pool[n].numParts; ++j) {
					int oi=pool[n].parts[j];
					if(oi!=i) {
						double dx=pop.getx(oi)-pop.getx(i);
						double dy=pop.gety(oi)-pop.gety(i);
						double dz=pop.getz(oi)-pop.getz(i);
						double dist=sqrt(dx*dx+dy*dy+dz*dz);
						double mag=pop.getTimeStep()*pop.getMass(oi)/(dist*dist*dist);
						pop.setvx(i,pop.getvx(i)+dx*mag);
						pop.setvy(i,pop.getvy(i)+dy*mag);
						pop.setvz(i,pop.getvz(i)+dz*mag);
					}
				}
			} else {
				double dx=pool[n].cmx-pop.getx(i);
				double dy=pool[n].cmy-pop.gety(i);
				double dz=pool[n].cmz-pop.getz(i);
				double dsqr=dx*dx+dy*dy+dz*dz;
				if(theta*theta*dsqr>pool[n].size*pool[n].size) {
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
					if(dsqr+pool[n].size*pool[n].size<minDist) return;
					doForce(pool[n].firstChild,pop,i);
					doForce(pool[n].firstChild+1,pop,i);
				}
			}
		}
		
		vector<KDNode> pool;
		int firstFree;
		int finalRoot;
		double theta;
		double minDist;
		int parallelRoot;
};

