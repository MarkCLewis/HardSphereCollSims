// FixedGridCollisionHash.h
// This class represents my standard method of creating a set of bins to
// quickly find pairs of particles that need to be tested for collisions.
// The grid size on it is a fixed grid where all bins are the same size.

#ifndef FIXED_GRID_COLLISION_HASH
#define FIXED_GRID_COLLISION_HASH

#ifdef PARALLEL
extern void printProcessNumber();
#endif

#include<vector>
#include<math.h>
#include "BinIndex.h"

class FixedGridCollisionHash {
	public:
		FixedGridCollisionHash() {
			buildsSinceLastAdjust=10000;
			first.resize(0);
			next.resize(0);
		}

		template<class Population>
		void build(Population &pop) {
			findBounds(pop);
			if(first.size()==0) {
				// This is the first time so build a hash in such a way as to
				// get roughly 10 particles per hash cell.  Then use the regular
				// algorithm to get a good size after that.
				double area=(maxx-minx)*(maxy-miny);
				minSpacing=sqrt(10.0*area/pop.getNumBodies());
				printf("area=%e, minSpacing=%e %e\n",area,minSpacing,10.0*area/pop.getNumBodies());
#ifndef PARALLEL
				buildReal(pop);
#endif
			}
#ifdef PARALLEL
			buildReal(pop);
#endif
			adjust(pop);
			buildReal(pop);
		}

		template<class Population>
		void buildReal(Population &pop) {
			//numBodies is used in parallelized for-loops to
			//conform to icc and OpenMP
			int numBodies = pop.getNumBodies();
#ifdef PARALLEL
			printProcessNumber();
#endif
			printf("bounds=%f %f %f %f\n",minx,maxx,miny,maxy);
			fflush(stdout);

			// At this point we should have a valid minSpacing value.
			// So we can use that to build the hash with.  Add 1 just
			// to be safe.
			factor=1.0/minSpacing;
			int numx=(int)((maxx-minx)*factor)+1;
			int numy=(int)((maxy-miny)*factor)+1;
			//fprintf(stderr, "%d gridSpacing = %e with Grid %d by %d\n", firstMade->getProcessNum(), minSpacing,numx, numy);
			//fflush(stderr);
#ifdef PARALLEL
			printProcessNumber();
#endif
			printf("Resize first to %d by %d\n",numx,numy);
			fflush(stdout);
			first.resize(numx);
			inUse.resize(numx);
			for(int i=0; i<numx; ++i) {
				first[i].resize(numy);
				inUse[i].resize(numy);
				for(int j=0; j<numy; ++j)
				{
				   first[i][j] = ParticleIndex{-1};
				   inUse[i][j] = false;
				}
			}
#ifdef PARALLEL
			printProcessNumber();
#endif
			printf("Hash resized to %d by %d for %d bodies\n",numx,numy,pop.getNumBodies());
			fflush(stdout);
			next.resize(pop.getNumBodies());
			xbin.resize(pop.getNumBodies());
			ybin.resize(pop.getNumBodies());
#ifdef PARALLEL
			printProcessNumber();
#endif
			printf("Done resizing\n");
			fflush(stdout);

			rowLock.resize(numx);
			for (int i = 0; i < (int)rowLock.size(); i++)
				omp_init_lock(&rowLock[i]);

			#pragma omp parallel for schedule(static)
			for(int i = 0; i<numBodies; ++i)
			{
				ParticleIndex pi = {i};
				BinIndex binx=calcBinX(pop,pi);
				BinIndex biny=calcBinY(pop,pi);
				if(binx.i<0) {
					binx.i=0;
				}
				if(binx.i>=(int)first.size()) {
					binx.i=first.size()-1;
				}
				if(biny.i<0) {
					biny.i=0;
				}
				if(biny.i>=(int)first[0].size()) {
					biny.i=first[0].size()-1;
				}
				xbin[pi.i]=binx;
				ybin[pi.i]=biny;

				omp_set_lock(&rowLock[binx.i]);
				next[i]=first[binx.i][biny.i];
				first[binx.i][biny.i]=pi;
				omp_unset_lock(&rowLock[binx.i]);
			}

			for (int i = 0; i < (int)rowLock.size(); i++)
				omp_destroy_lock(&rowLock[i]);

			++buildsSinceLastAdjust;
		}

		int getMaxX() {
			return first.size();
		}

		int getMaxY() {
			return first[0].size();
		}

		double getMinSpacing() {
			return 3*minSpacing;
		}

		ParticleIndex getFirst(BinIndex x,BinIndex y) {
			return first[x.i][y.i];
		}

		ParticleIndex getNext(ParticleIndex part) {
			return next[part.i];
		}

		template<class Population>
		BinIndex calcBinX(const Population &pop,ParticleIndex pi) {
			return BinIndex{(int)((pop.getx(pi)-minx)*factor)};
		}

		template<class Population>
		BinIndex calcBinY(Population &pop,ParticleIndex pi) {
			return BinIndex{(int)((pop.gety(pi)-miny)*factor)};
		}

		BinIndex getBinX(ParticleIndex pi) {
			return xbin[pi.i];
		}

		BinIndex getBinY(ParticleIndex pi) {
			return ybin[pi.i];
		}

		bool isSafe(ParticleIndex pi)
		{
			BinIndex x=xbin[pi.i];
			BinIndex y=ybin[pi.i];
			int startx = (x.i>0) ? -1 : 0;
			int endx = (x.i<(int)inUse.size()-1) ? 2 : 1;
			int starty = (y.i>0) ? -1 : 0;
			int endy = (y.i<(int)inUse[0].size()-1) ? 2 : 1;
			for (int i = startx; i < endx; i++)
			{
				for (int j = starty; j < endy; j++)
				{
					if (inUse[x.i+i][y.i+j])
						return false;
					/*else {
						printf("%d\t (%d,%d) is safe\n",omp_get_thread_num(),x+i,y+j);
						fflush(stdout);
					}*/
				}
			}
			return true;
		}

		void setInUse(ParticleIndex pi, bool val)
		{
			inUse[xbin[pi.i].i][ybin[pi.i].i] = val;
		}

	private:
		template<class Population>
		void findBounds(Population &pop) {
			int numBodies = pop.getNumBodies();
			minx=pop.getx(ParticleIndex{0});
			maxx=pop.getx(ParticleIndex{0});
			miny=pop.gety(ParticleIndex{0});
			maxy=pop.gety(ParticleIndex{0});

			#pragma omp parallel
			{
				double lminx=minx, lmaxx=maxx;
				double lminy=miny, lmaxy=maxy;
				#pragma omp for schedule(static)
				for(int i=1; i<numBodies; ++i)
				{
					ParticleIndex pi = {i};
					if (pop.getx(pi)<lminx) lminx=pop.getx(pi);
					if (pop.getx(pi)>lmaxx) lmaxx=pop.getx(pi);
					if (pop.gety(pi)<lminy) lminy=pop.gety(pi);
					if (pop.gety(pi)>lmaxy) lmaxy=pop.gety(pi);
				}
				#pragma omp critical
				{
					if (lminx<minx) minx=lminx;
					if (lmaxx>maxx) maxx=lmaxx;
					if (lminy<miny) miny=lminy;
					if (lmaxy>maxy) maxy=lmaxy;
				}
			}
		}
		// This function assigns a safe value to minSpacing for the
		// hash to use.
		template<class Population>
		void adjust(Population &pop) {
			buildsSinceLastAdjust=0;
			if(pop.getNumBodies()<1000) return;
			// Run through calculating a velocity dispersion.  I normally
			// use just RMS, but it might be safer to use something that
			// is more sensative to the higher relative velocities.
			double sum=0.0,cnt=0.0;
			double max=pop.getNumBodies()/100;
			while(cnt<max) {
				int x=lrand48()%first.size();
				int y=lrand48()%first[0].size();
				ParticleIndex f=first[x][y];
//				printf("%d %d %d %d %d %d\n",f,cnt,x,y,first.size(),first[0].size());
				if(f.i>-1) {
					ParticleIndex n=next[f.i];
					while(n.i>-1) {
						// Increment sum.
						double dvx=pop.getvx(f)-pop.getvx(n);
						double dvy=pop.getvy(f)-pop.getvy(n);
						double dvz=pop.getvz(f)-pop.getvz(n);
						sum+=dvx*dvx+dvy*dvy+dvz*dvz;
						cnt+=1.0;
						n=next[n.i];
					}
					int x2,y2;
					n=f;
					while(n.i==f.i) {
						x2=x+(lrand48()%3)-1;
						y2=y+(lrand48()%3)-1;
						if(x2>=0 && x2<(int)first.size() && y2>=0 && y2<(int)first[x2].size())
							n=first[x2][y2];
						else n.i=-1;
					}
					while(n.i>-1) {
						// Increment sum.
						double dvx=pop.getvx(f)-pop.getvx(n);
						double dvy=pop.getvy(f)-pop.getvy(n);
						double dvz=pop.getvz(f)-pop.getvz(n);
						sum+=dvx*dvx+dvy*dvy+dvz*dvz;
						cnt+=1.0;
						n=next[n.i];
					}
				}
			}
			minSpacing=2.0*pop.getMaxParticleRadius()+5.0*sqrt(sum/cnt)*pop.getTimeStep();
			printf("minSpacing1=%e\n",minSpacing);
			double area=(maxx-minx)*(maxy-miny);
			if(area/(minSpacing*minSpacing)>maxRatio*pop.getNumBodies()) {
				minSpacing=sqrt(area/(maxRatio*pop.getNumBodies()));
			}
			printf("minSpacing2=%e\n",minSpacing);
			//fprintf(stderr, "%d cnt=%e sum=%e\n", firstMade->getProcessNum(), cnt, sum);
			//fflush(stderr);
		}

		std::vector<std::vector<ParticleIndex> > first;
		std::vector<ParticleIndex> next;
		std::vector<BinIndex> xbin;
		std::vector<BinIndex> ybin;
		int buildsSinceLastAdjust;
		double minSpacing;
		double factor;
		double minx,maxx,miny,maxy;
		const static double maxRatio;
		std::vector<std::vector<bool> > inUse;
		std::vector<omp_lock_t> rowLock;
};

const double FixedGridCollisionHash::maxRatio=7.0;

#endif
