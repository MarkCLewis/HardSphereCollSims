#ifndef PARTICLE_2D_BINNED_OUTPUT
#define PARTICLE_2D_BINNED_OUTPUT

#include <cstdio>
#include <vector>
#include <algorithm>

using std::vector;
using std::fwrite;
using std::sort;

struct PartBin {
	PartBin():count(0),x(0),y(0),minx(0),maxx(0),tau(0),e(0),i(0),phi(0) {}
	void zero() { count=0; x=0; y=0; minx=0; maxx=0; tau=0; e=0; i=0; phi=0; }
	int count;
	float x;
	float y;
	float minx;
	float maxx;
	float tau;
	float e;
	float i;
	float phi;
};

template<class GCType>
struct CoordsR {
	CoordsR() {}
	CoordsR(CartCoords nc,GCType ngc,double nr):c(nc),gc(ngc),r(nr) {}
	CartCoords c;
	GCType gc;
	double r;
};

template<class GCType>
bool cartCompFunc(const CoordsR<GCType> &p1,const CoordsR<GCType> &p2) {
	if(p1.c.x<p2.c.x) return true;
	return false;
}

template<class GCType>
bool gcCompFunc(const CoordsR<GCType> &p1,const CoordsR<GCType> &p2) {
	if(p1.gc.X<p2.gc.X) return true;
	return false;
}

template<class BoundaryCondition,class GCType>
class SpecificParticle2DBinnedOutput {
	public:
		SpecificParticle2DBinnedOutput(int oi,int bcx,int bcy,BoundaryCondition &bc):interval(oi),
				binCountX(bcx),binCountY(bcy),bounds(bc) {
			cartParts.resize(binCountY);
			gcParts.resize(binCountY);
			cartBins.resize(binCountY);
			gcBins.resize(binCountY);
			for(int i=0; i<binCountY; ++i) {
				cartBins[i].resize(binCountX);
				gcBins[i].resize(binCountX);
			}
			firstFree=-1;
		}

		bool doOnStep(int step) {
			return step%interval==0;
		}

		void init(int step) {
			miny=bounds.getMinYTotal();
			maxy=bounds.getMaxYTotal();

			for(int i=0; i<binCountY; ++i) {
				cartParts[i]=-1;
				gcParts[i]=-1;
				for(int j=0; j<binCountX; ++j) {
					cartBins[i][j].count=0;
					cartBins[i][j].x=0;
					cartBins[i][j].y=0;
					cartBins[i][j].minx=0;
					cartBins[i][j].maxx=0;
					cartBins[i][j].tau=0;
					cartBins[i][j].e=0;
					cartBins[i][j].i=0;
					cartBins[i][j].phi=0;
					gcBins[i][j].count=0;
					gcBins[i][j].x=0;
					gcBins[i][j].y=0;
					gcBins[i][j].minx=0;
					gcBins[i][j].maxx=0;
					gcBins[i][j].tau=0;
					gcBins[i][j].e=0;
					gcBins[i][j].i=0;
					gcBins[i][j].phi=0;
				}
			}
			vector<bool> isFree(pool.size());
			for(int i=0; i<isFree.size(); ++i) {
				isFree[i]=false;
			}
			int n=firstFree;
			while(n>-1) {
				isFree[n]=true;
				n=(int)pool[n][0].r;
			}
			for(int i=0; i<isFree.size(); ++i) {
				if(!isFree[i]) {
					printf("LEAK!!! %d\n",i,pool[i].size());
				}
			}
		}

		void process(int num,vector<CartCoords> &c,vector<double> &r) {
			printf("Process a machine.\n");
			minCartBin=1000000;
			maxCartBin=0;
			minGCBin=1000000;
			maxGCBin=0;
			for(int i=0; i<num; ++i) {
				GCType gc(c[i]);
				int cartBin=(int)((c[i].y-miny)*binCountY/(maxy-miny));
				int gcBin=(int)((gc.Y-miny)*binCountY/(maxy-miny));
				if(cartBin<0) cartBin=0;
					else if(cartBin>=binCountY) cartBin=binCountY-1;
				if(gcBin<0) gcBin=0;
					else if(gcBin>=binCountY) gcBin=binCountY-1;
				if(cartBin<minCartBin) minCartBin=cartBin;
				if(cartBin>maxCartBin) maxCartBin=cartBin;
				if(gcBin<minGCBin) minGCBin=gcBin;
				if(gcBin>maxGCBin) maxGCBin=gcBin;
				CoordsR<GCType> cr(c[i],gc,r[i]);
				if(cartParts[cartBin]==-1) {
					cartParts[cartBin]=grabPoolRow();
				}
				pool[cartParts[cartBin]].push_back(cr);
				if(gcParts[gcBin]==-1) {
					gcParts[gcBin]=grabPoolRow();
				}
				pool[gcParts[gcBin]].push_back(cr);
			}

			// Process and clear below min.
			for(int i=0; i<minCartBin; ++i) {
				if(cartParts[i]>=0) {
					sortAndBinCart(i);
				}
			}
			for(int i=0; i<minGCBin; ++i) {
				if(gcParts[i]>=0) {
					sortAndBinGC(i);
				}
			}
			//printf("Done with processing\n");
			//fflush(stdout);
		}

		int grabPoolRow(void) {
			if(firstFree==-1) {
				pool.resize(pool.size()+1);
				//printf("Made new pool col %d\n",pool.size());
				return pool.size()-1;
			} else {
				int tmp=firstFree;
				//printf("Grabbing free %d\n",firstFree);
				firstFree=(int)pool[firstFree][0].r;
				pool[tmp].clear();
				return tmp;
			}
		}

		void sortAndBinCart(int i) {
			sort(pool[cartParts[i]].begin(),pool[cartParts[i]].end(),cartCompFunc<GCType>);
			for(int j=0; j<binCountX; ++j) {
				PartBin &bin=cartBins[i][j];
				unsigned int start=j*pool[cartParts[i]].size()/binCountX;
				unsigned int end=(j+1)*pool[cartParts[i]].size()/binCountX;
				bin.minx=pool[cartParts[i]][start].c.x;
				bin.maxx=pool[cartParts[i]][start].c.x;
				for(unsigned int k=start; k<end && k<pool[cartParts[i]].size(); ++k) {
					bin.count++;
					double x=pool[cartParts[i]][k].c.x;
					bin.x+=x;
					bin.y+=pool[cartParts[i]][k].c.y;
					if(x<bin.minx) bin.minx=x;
					if(x>bin.maxx) bin.maxx=x;
					bin.tau+=3.14159*pool[cartParts[i]][k].r*pool[cartParts[i]][k].r;
					bin.e+=pool[cartParts[i]][k].gc.e;
					bin.i+=pool[cartParts[i]][k].gc.i;
					float phi=pool[cartParts[i]][k].gc.phi;
					if(bin.count>1) {
						while(bin.phi/(bin.count-1)-phi>3.14159) {
							phi+=2.0*3.14159;
						} 
						while(bin.phi/(bin.count-1)-phi<-3.14159) {
							phi-=2.0*3.14159;
						}
					}
					bin.phi+=phi;
				}
				if(bin.count>0) {
					bin.x/=bin.count;
					bin.y/=bin.count;
					double xmin=bin.minx,xmax=bin.maxx;
					if(start>0) xmin=0.5*(xmin+pool[cartParts[i]][start-1].c.x);
					if(end<pool[gcParts[i]].size()) xmax=0.5*(xmax+pool[cartParts[i]][end].c.x);
					bin.tau*=binCountY/((xmax-xmin)*(maxy-miny));
					bin.e/=bin.count;
					bin.i/=bin.count;
					bin.phi/=bin.count;
				} else {
					bin.x=pool[cartParts[i]][start].c.x;
					bin.y=pool[cartParts[i]][start].c.y;
				}
			}
			pool[cartParts[i]].resize(1);
			pool[cartParts[i]][0].r=firstFree;
			firstFree=cartParts[i];
			cartParts[i]=-1;
			//printf("Cart pool freed %d\n",firstFree);
		}

		void sortAndBinGC(int i) {
			//printf("Sort and bin gc for %d %d %d\n",i,gcParts[i],pool[gcParts[i]].size());
			sort(pool[gcParts[i]].begin(),pool[gcParts[i]].end(),gcCompFunc<GCType>);
			for(int j=0; j<binCountX; ++j) {
				PartBin &bin=gcBins[i][j];
				unsigned int start=j*pool[gcParts[i]].size()/binCountX;
				unsigned int end=(j+1)*pool[gcParts[i]].size()/binCountX;
				bin.minx=pool[gcParts[i]][start].gc.X;
				bin.maxx=pool[gcParts[i]][start].gc.X;
				for(unsigned int k=start; k<end && k<pool[gcParts[i]].size(); ++k) {
					bin.count++;
					double x=pool[gcParts[i]][k].gc.X;
					bin.x+=x;
					bin.y+=pool[gcParts[i]][k].gc.Y;
					if(x<bin.minx) bin.minx=x;
					if(x>bin.maxx) bin.maxx=x;
					bin.tau+=3.14159*pool[gcParts[i]][k].r*pool[gcParts[i]][k].r;
					bin.e+=pool[gcParts[i]][k].gc.e;
					bin.i+=pool[gcParts[i]][k].gc.i;
					float phi=pool[gcParts[i]][k].gc.phi;
					if(bin.count>1) {
						while(bin.phi/(bin.count-1)-phi>3.14159) {
							phi+=2.0*3.14159;
						} 
						while(bin.phi/(bin.count-1)-phi<-3.14159) {
							phi-=2.0*3.14159;
						}
					}
					bin.phi+=phi;
				}
				if(bin.count>0) {
					bin.x/=bin.count;
					bin.y/=bin.count;
					double xmin=bin.minx,xmax=bin.maxx;
					if(start>0) xmin=0.5*(xmin+pool[gcParts[i]][start-1].gc.X);
					if(end<pool[gcParts[i]].size()) xmax=0.5*(xmax+pool[gcParts[i]][end].gc.X);
					bin.tau*=binCountY/((xmax-xmin)*(maxy-miny));
					bin.e/=bin.count;
					bin.i/=bin.count;
					bin.phi/=bin.count;
				} else {
					bin.x=pool[gcParts[i]][start].gc.X;
					bin.y=pool[gcParts[i]][start].gc.Y;
				}
			}
			pool[gcParts[i]].resize(1);
			pool[gcParts[i]][0].r=firstFree;
			firstFree=gcParts[i];
			gcParts[i]=-1;
			//printf("GC pool freed %d\n",firstFree);
		}

		void finalize(int step) {
			printf("Finalize %d\n",pool.size());
			// Process and clear remaining.
			for(int i=0; i<cartParts.size(); ++i) {
				if(cartParts[i]>=0) {
					sortAndBinCart(i);
				}
			}
			for(int i=0; i<gcParts.size(); ++i) {
				if(gcParts[i]>=0) {
					sortAndBinGC(i);
				}
			}

			char buf[100];
			sprintf(buf,"ParticleCartBinned.%d.bin",step);
			FILE *cartFout=writeHeader(buf);
			sprintf(buf,"ParticleGCBinned.%d.bin",step);
			FILE *gcFout=writeHeader(buf);
			for(int i=0; i<binCountY; ++i) {
				writeLineHeader(cartFout,miny+(i+0.5)*(maxy-miny)/binCountY);
				writeLineHeader(gcFout,miny+(i+0.5)*(maxy-miny)/binCountY);
				for(int j=0; j<binCountX; ++j) {
					writeBin(cartFout,cartBins[i][j]);
				}
				for(int j=0; j<binCountX; ++j) {
					writeBin(gcFout,gcBins[i][j]);
				}
			}
			fclose(cartFout);
			fclose(gcFout);
			cartFout=0;
			gcFout=0;
		}
	private:
		FILE *writeHeader(char *filename) {
			FILE *fout=fopen(filename,"wb");
		
			// Write out header information.
			fwrite(&binCountX,sizeof(int),1,fout);
			fwrite(&binCountY,sizeof(int),1,fout);
			int numParams=1,numValues=8;
			fwrite(&numParams,sizeof(int),1,fout);
			fwrite("Count",sizeof(char),6,fout);
			fwrite(&numValues,sizeof(int),1,fout);
			fwrite("x",sizeof(char),2,fout);
			fwrite("y",sizeof(char),2,fout);
			fwrite("minx",sizeof(char),5,fout);
			fwrite("maxx",sizeof(char),5,fout);
			fwrite("tau",sizeof(char),4,fout);
			fwrite("e",sizeof(char),2,fout);
			fwrite("i",sizeof(char),2,fout);
			fwrite("phi",sizeof(char),4,fout);
			return fout;
		}
		
		void writeLineHeader(FILE* fout,float y) {
			fwrite(&y,sizeof(float),1,fout);
		}
		
		void writeBin(FILE *fout,PartBin &data) {
			fwrite(&(data.count),sizeof(int),1,fout);
			fwrite(&(data.x),sizeof(float),1,fout);
			fwrite(&(data.y),sizeof(float),1,fout);
			fwrite(&(data.minx),sizeof(float),1,fout);
			fwrite(&(data.maxx),sizeof(float),1,fout);
			fwrite(&(data.tau),sizeof(float),1,fout);
			fwrite(&(data.e),sizeof(float),1,fout);
			fwrite(&(data.i),sizeof(float),1,fout);
			fwrite(&(data.phi),sizeof(float),1,fout);
		}
		
		int interval;
		int binCountX;
		int binCountY;
		BoundaryCondition &bounds;
	
		vector<vector<CoordsR<GCType> > > pool;
		int firstFree;
		vector<int> cartParts;
		vector<int> gcParts;
		vector<vector<PartBin> > cartBins;
		vector<vector<PartBin> > gcBins;
		int minCartBin,maxCartBin;
		int minGCBin,maxGCBin;
		double miny,maxy;
		
};

#endif
