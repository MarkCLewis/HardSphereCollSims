/**
 * This header file includes the code to do a forcing based on a tree.
**/

#include <vector>
#include <algorithm>

using std::vector;
using std::min;
using std::max;

const int MAX_PARTS=2;

struct QuadNode {
	double cx,cy,size;
	int firstChild;
	int parts[MAX_PARTS];
	int numParts;
	double cmx,cmy,cmz,mass;
	
	int childNum(double x,double y) {
		int ret=0;
		if(x>cx) ret|=1;
		if(y>cy) ret|=2;
		return ret;
	}
};

class QuadGravTree {
	public:
		QuadGravTree(double t):theta(t) {}
	
		template<class Population>
		void applyForce(Population &pop) {
			printf("Do tree gravity.\n");
			pool.resize(4*pop.getNumBodies());
			pool[0].numParts=0;
			firstFree=1;
			double minx=pop.getx(0),maxx=pop.getx(0);
			double miny=pop.gety(0),maxy=pop.gety(0);
			for(int i=1; i<pop.getNumBodies(); ++i) {
				if(minx>pop.getx(i)) minx=pop.getx(i);
				if(maxx<pop.getx(i)) maxx=pop.getx(i);
				if(miny>pop.gety(i)) miny=pop.gety(i);
				if(maxy<pop.gety(i)) maxy=pop.gety(i);
			}
			pool[0].size=max(maxx-minx,maxy-miny);
			pool[0].cx=0.5*(maxx+minx);
			pool[0].cy=0.5*(maxy+miny);
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
				addParticle(pool[n].firstChild+pool[n].childNum(pop.getx(i),pop.gety(i)),pop,i);
			} else if(pool[n].numParts<MAX_PARTS) {
				pool[n].parts[pool[n].numParts]=i;
				++pool[n].numParts;
			} else {
				pool[n].firstChild=firstFree;
				pool[firstFree].numParts=0;
				pool[firstFree+1].numParts=0;
				pool[firstFree+2].numParts=0;
				pool[firstFree+3].numParts=0;
				double hsize=pool[n].size*0.5;
				double qsize=pool[n].size*0.25;
				pool[firstFree].size=hsize;
				pool[firstFree+1].size=hsize;
				pool[firstFree+2].size=hsize;
				pool[firstFree+3].size=hsize;
				pool[firstFree].cx=pool[n].cx-qsize;
				pool[firstFree+1].cx=pool[n].cx+qsize;
				pool[firstFree+2].cx=pool[n].cx-qsize;
				pool[firstFree+3].cx=pool[n].cx+qsize;
				pool[firstFree].cy=pool[n].cy-qsize;
				pool[firstFree+1].cy=pool[n].cy-qsize;
				pool[firstFree+2].cy=pool[n].cy+qsize;
				pool[firstFree+3].cy=pool[n].cy+qsize;
				firstFree+=4;
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
				for(int j=0; j<4; ++j) {
					finalize(pool[n].firstChild+j,pop);
					pool[n].mass+=pool[pool[n].firstChild+j].mass;
					pool[n].cmx+=pool[pool[n].firstChild+j].mass*pool[pool[n].firstChild+j].cmx;
					pool[n].cmy+=pool[pool[n].firstChild+j].mass*pool[pool[n].firstChild+j].cmy;
					pool[n].cmz+=pool[pool[n].firstChild+j].mass*pool[pool[n].firstChild+j].cmz;
				}
			} else {
				for(int j=0; j<pool[n].numParts; ++j) {
					pool[n].mass+=pop.getMass(pool[n].parts[j]);
					pool[n].cmx+=pop.getMass(pool[n].parts[j])*pop.getx(pool[n].parts[j]);
					pool[n].cmy+=pop.getMass(pool[n].parts[j])*pop.gety(pool[n].parts[j]);
					pool[n].cmz+=pop.getMass(pool[n].parts[j])*pop.getz(pool[n].parts[j]);
				}
			}
			if(pool[n].mass>0.0) {
				pool[n].cmx/=pool[n].mass;
				pool[n].cmy/=pool[n].mass;
				pool[n].cmz/=pool[n].mass;
			}
		}

		template<class Population>
		void doForce(int n,Population &pop,int i) {
			if(pool[n].mass<=0.0) return;
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
					doForce(pool[n].firstChild,pop,i);
					doForce(pool[n].firstChild+1,pop,i);
					doForce(pool[n].firstChild+2,pop,i);
					doForce(pool[n].firstChild+3,pop,i);
				}
			}
		}
		
		vector<QuadNode> pool;
		int firstFree;
		double theta;
};

struct KDNode {
	int splitDim;
	double splitVal;
	int firstChild;
	int parts[MAX_PARTS];
	int numParts;
	
	double cmx,cmy,cmz,mass;
	double minx,maxx,miny,maxy,size;
	
	int childNum(double x,double y) {
		int ret=0;
		if((splitDim==0 && x>splitVal) || (splitDim==1 && y>splitVal)) ret=1;
		return ret;
	}
};

struct ForceInfo {
	int n,part;
};

class KDGravTree {
	public:
		KDGravTree(double t):theta(t) {
			maxDepth = (int) (log(omp_get_max_threads())/log(2) + 2);
		}
	
		template<class Population>
		void applyForce(Population &pop) {
			printf("Do tree gravity on %d bodies.\n",pop.getNumBodies());
			pool.resize(2*pop.getNumBodies());
			pool[0].numParts=0;
			firstFree=1;
			printf("Build tree on %d bodies.\n",pop.getNumBodies());	
			for(int i=0; i<pop.getNumBodies(); ++i) {
				addParticle(0,pop,i);
			}

			printf("Finalize tree on %d bodies.\n",pop.getNumBodies());
			finalize(0,pop);
			printf("Do forcing on %d bodies\n",pop.getNumBodies());
			fflush(stdout);

			double start = omp_get_wtime();
			int nb = pop.getNumBodies();
			#pragma omp parallel for schedule(dynamic,100)
			for (int i=0; i<nb; i++) {
				doForce(0,pop,i);
				pop.adjustAfterForce(i);
			}
			fprintf(stderr,"Forcing took %g seconds.\n",omp_get_wtime()-start);
			fflush(stderr);
			printf("Done with tree force on %d bodies.\n",pop.getNumBodies());
			fflush(stdout);
		}
	private:
		template<class Population>
		void addParticle(int n,Population &pop,int i) {
			if(pool[n].numParts<0) {
				addParticle(pool[n].firstChild+pool[n].childNum(pop.getx(i),pop.gety(i)),pop,i);
			} else if(pool[n].numParts<MAX_PARTS) {
				pool[n].parts[pool[n].numParts]=i;
				++pool[n].numParts;
			} else {
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
				for(int j=0; j<2; ++j) {
					finalize(pool[n].firstChild+j,pop);
					pool[n].mass+=pool[pool[n].firstChild+j].mass;
					pool[n].cmx+=pool[pool[n].firstChild+j].mass*pool[pool[n].firstChild+j].cmx;
					pool[n].cmy+=pool[pool[n].firstChild+j].mass*pool[pool[n].firstChild+j].cmy;
					pool[n].cmz+=pool[pool[n].firstChild+j].mass*pool[pool[n].firstChild+j].cmz;
				}
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
		void doForceTop(int n,Population &pop,int i,vector<ForceInfo> &info, int depth=0) {
			if (depth >= maxDepth) {
				info.resize(info.size()+1);
				info[info.size()-1].n = n;
				info[info.size()-1].part = i;
				return;
			}
			if(pool[n].mass<=0.0) return;
			if(n<0 || n>=(int)pool.size()) printf("BAD! %d\n",n);
			if(i<0 || i>=pop.getNumBodies()) printf("BAD2! %d\n",i);
			if(pool[n].numParts>0) {
				for(int j=0; j<pool[n].numParts; ++j) {
					int oi=pool[n].parts[j];
//					if(pop.getNumBodies()==411581) {
//						printf("oi=%d\n",oi);
//						fflush(stdout);
//					}
					if(oi<0 || oi>=pop.getNumBodies()) printf("BAD3! %d\n",oi);
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
//					if(pop.getNumBodies()==411581) {
//						printf("interact CM\n");
//						fflush(stdout);
//					}
					double dist=sqrt(dsqr);
					double mag=pop.getTimeStep()*pool[n].mass/(dist*dist*dist);
					pop.setvx(i,pop.getvx(i)+dx*mag);
					pop.setvy(i,pop.getvy(i)+dy*mag);
					pop.setvz(i,pop.getvz(i)+dz*mag);
				} else {
//					if(pop.getNumBodies()==411581) {
//						printf("recurse 1\n");
//						fflush(stdout);
//					}
					doForceTop(pool[n].firstChild,pop,i,info,depth+1);
//					if(pop.getNumBodies()==411581) {
//						printf("recurse 2\n");
//						fflush(stdout);
//					}
					doForceTop(pool[n].firstChild+1,pop,i,info,depth+1);
				}
			}
//			if(pop.getNumBodies()==411581) {
//				printf("exit force %d %d\n",n,i);
//				fflush(stdout);
//			}
		}

		template <class Population>
		void doForce(int n,Population &pop,int i) {
//			if(pop.getNumBodies()==411581) {
//				printf("%d %d\n",n,i);
//				fflush(stdout);
//			}
			if(pool[n].mass<=0.0) return;
			if(n<0 || n>=(int)pool.size()) printf("BAD! %d\n",n);
			if(i<0 || i>=pop.getNumBodies()) printf("BAD2! %d\n",i);
			if(pool[n].numParts>0) {
				for(int j=0; j<pool[n].numParts; ++j) {
					int oi=pool[n].parts[j];
//					if(pop.getNumBodies()==411581) {
//						printf("oi=%d\n",oi);
//						fflush(stdout);
//					}
					if(oi<0 || oi>=pop.getNumBodies()) printf("BAD3! %d\n",oi);
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
//					if(pop.getNumBodies()==411581) {
//						printf("interact CM\n");
//						fflush(stdout);
//					}
					double dist=sqrt(dsqr);
					double mag=pop.getTimeStep()*pool[n].mass/(dist*dist*dist);
					pop.setvx(i,pop.getvx(i)+dx*mag);
					pop.setvy(i,pop.getvy(i)+dy*mag);
					pop.setvz(i,pop.getvz(i)+dz*mag);
				} else {
//					if(pop.getNumBodies()==411581) {
//						printf("recurse 1\n");
//						fflush(stdout);
//					}
					doForce(pool[n].firstChild,pop,i);
//					if(pop.getNumBodies()==411581) {
//						printf("recurse 2\n");
//						fflush(stdout);
//					}
					doForce(pool[n].firstChild+1,pop,i);
				}
			}
//			if(pop.getNumBodies()==411581) {
//				printf("exit force %d %d\n",n,i);
//				fflush(stdout);
//			}
		}

		
		vector<KDNode> pool;
		int firstFree;
		double theta;
		int maxDepth;
};

