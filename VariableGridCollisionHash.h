// VariableGridCollisionHash.h
// This class represents a method of binning up particles to search for
// collisions that is adaptive yet still provides O(1) placement and retrievals.
// The idea is that there is a simple mapping for a very high resolution 1-D
// array in both x and y to the cells of the more course grained 2-D grid that
// is used to search for collisions.  This allows the total number of grid cells
// to stay

#ifndef VARIABLE_GRID_COLLISION_HASH
#define VARIABLE_GRID_COLLISION_HASH

#include<vector>
#include<math.h>

class VariableGridCollisionHash {
	public:
		VariableGridCollisionHash() {
			buildsSinceLastAdjust=10000;
			first.resize(0);
			next.resize(0);
		}

		template<class Population>
		void build(Population &pop) {

			minx=pop.getx(0);
			maxx=pop.getx(0);
			miny=pop.gety(0);
			maxy=pop.gety(0);
			for(int i=1; i<pop.getNumBodies(); ++i) {
				if(pop.getx(i)<minx) minx=pop.getx(i);
				if(pop.getx(i)>maxx) maxx=pop.getx(i);
				if(pop.gety(i)<miny) miny=pop.gety(i);
				if(pop.gety(i)>maxy) maxy=pop.gety(i);
			}

			if(buildsSinceLastAdjust>1000) {
				adjust(pop);
			}
			// At this point we should have a valid minSpacing and factor value.
			// So we can use that to build the hash with.  Add 1 just
			// to be safe.  First resize all the vector.
			int numx=(int)((maxx-minx)*factor)+1;
			int numy=(int)((maxy-miny)*factor)+1;
			xMap.resize(numx);
			for(int i=0; i<numx; ++i) xMap[i]=0;
			yMap.resize(numy);
			for(int i=0; i<numy; ++i) yMap[i]=0;
			int sy=(int)(sqrt(maxRatio*pop.getNumBodies()*(maxy-miny)/(maxx-minx)));
			int sx=(int)(maxRatio*pop.getNumBodies()/sy);
			first.resize(sx+1);
			for(unsigned int i=0; i<first.size(); ++i) {
				first[i].resize(sy+1);
				for(unsigned int j=0; j<first[i].size(); ++j) first[i][j]=-1;
			}
			next.resize(pop.getNumBodies());

			// Build mappings.
			xBin.resize(pop.getNumBodies());
			yBin.resize(pop.getNumBodies());
			for(int i=0; i<pop.getNumBodies(); ++i) {
				int binx=(int)((pop.getx(i)-minx)*factor);
				xMap[binx]++;
				xBin[i]=binx;
				int biny=(int)((pop.gety(i)-miny)*factor);
				yMap[biny]++;
				yBin[i]=biny;
//				if(i<100) printf("%d %e   %e %e->%d  %e %e->%d\n",i,minSpacing,pop.getx(i),minx,binx,pop.gety(i),miny,biny);
			}

			int pTot=pop.getNumBodies()/sx;
			int pCnt=0;
			int binCnt=0;
			for(int i=0; i<numx; ++i) {
				pCnt+=xMap[i];
				xMap[i]=binCnt;
//				printf("xmap %d=%d  %d %d\n",i,binCnt,pCnt,pTot);
				if(pCnt>=pTot) {
					binCnt++;
					pCnt=0;
				}
			}
			pTot=pop.getNumBodies()/sy;
			pCnt=0;
			binCnt=0;
			for(int i=0; i<numy; ++i) {
				pCnt+=yMap[i];
				yMap[i]=binCnt;
//				printf("ymap %d=%d  %d %d\n",i,binCnt,pCnt,pTot);
				if(pCnt>=pTot) {
					binCnt++;
					pCnt=0;
				}
			}

			// Now build the actual hash with the variable sizes given by the
			// maps.
			for(int i=0; i<pop.getNumBodies(); ++i) {
				int binx=xBin[i];
				int biny=yBin[i];
				next[i]=first[xMap[binx]][yMap[biny]];
				first[xMap[binx]][yMap[biny]]=i;
			}

			++buildsSinceLastAdjust;
		}

		int getMaxX() {
			return first.size();
		}

		int getMaxY() {
			return first[0].size();
		}

		double getMinSpacing() {
			return minSpacing;
		}

		int getFirst(int x,int y) {
			return first[x][y];
		}

		int getNext(int part) {
			return next[part];
		}

		template<class Population>
		int getBinX(const Population &pop,int i) {
			int binx=(int)((pop.getx(i)-minx)*factor);
			return xMap[binx];
		}

		template<class Population>
		int getBinY(Population &pop,int i) {
			int biny=(int)((pop.gety(i)-miny)*factor);
			return yMap[biny];
		}
	private:
		// This function assigns a safe value to minSpacing for the
		// hash to use.
		template<class Population>
		void adjust(Population &pop) {
			buildsSinceLastAdjust=0;
			if(first.size()==0) {
				// This is the first time so build a hash in such a way as to
				// get roughly 10 particles per hash cell.  Then use the regular
				// algorithm to get a good size after that.
				double area=(maxx-minx)*(maxy-miny);
				minSpacing=sqrt(10.0*area/pop.getNumBodies());
				factor=1.0/minSpacing;
				build(pop);
			}
			// Run through calculating a velocity dispersion.  I normally
			// use just RMS, but it might be safer to use something that
			// is more sensative to the higher relative velocities.
			double sum=0.0,cnt=0.0;
			while(cnt<pop.getNumBodies()/100.0) {
				int x=lrand48()%first.size();
				int y=lrand48()%first[0].size();
//				printf("%d %d %d %d\n",x,y,first.size(),first[0].size());
				int f=first[x][y];
				if(f>-1) {
//					printf("%d %d %d %d\n",x,y,f,next.size());
					int n=next[f];
					while(n>-1) {
						// Increment sum.
						double dvx=pop.getvx(f)-pop.getvx(n);
						double dvy=pop.getvy(f)-pop.getvy(n);
						double dvz=pop.getvz(f)-pop.getvz(n);
						sum+=dvx*dvx+dvy*dvy+dvz*dvz;
						cnt+=1.0;
						n=next[n];
					}
					int x2,y2;
					n=f;
					while(n==f) {
						x2=x+(lrand48()%3)-1;
						y2=y+(lrand48()%3)-1;
						if(x2>0 && x2<(int)first.size() && y2>0 && y2<(int)first.size())
							n=first[x2][y2];
						else n=-1;
					}
					while(n>-1) {
						// Increment sum.
						double dvx=pop.getvx(f)-pop.getvx(n);
						double dvy=pop.getvy(f)-pop.getvy(n);
						double dvz=pop.getvz(f)-pop.getvz(n);
						sum+=dvx*dvx+dvy*dvy+dvz*dvz;
						cnt+=1.0;
						n=next[n];
					}
				}
			}
			minSpacing=2.0*pop.getMaxParticleRadius()+3.0*sqrt(sum/cnt)*pop.getTimeStep();
			factor=1.0/minSpacing;
		}

		std::vector<int> xMap;
		std::vector<int> yMap;
		std::vector<int> xBin;
		std::vector<int> yBin;
		std::vector<std::vector<int> > first;
		std::vector<int> next;
		int buildsSinceLastAdjust;
		double minSpacing;
		double factor;
		double minx,maxx,miny,maxy;
		const static double maxRatio;
};

const double VariableGridCollisionHash::maxRatio=7.0;

#endif
