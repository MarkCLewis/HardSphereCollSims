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

template<class Population>
class CartCompareFunc {
	public:
		CartCompareFunc(Population &p):pop(p) {}
		bool operator()(int p1,int p2) {
			double x1=pop.getx(p1);
			double x2=pop.getx(p2);
			if(x1<x2) return true;
			return false;
		}
	private:
		Population &pop;
};

template<class Population>
class GCCompareFunc {
	public:
		GCCompareFunc(Population &p):pop(p) {}
		bool operator()(int p1,int p2) {
			double x1=pop.getX(p1);
			double x2=pop.getX(p2);
			if(x1<x2) return true;
			return false;
		}
	private:
		Population &pop;
};

class Particle2DBinnedOutput {
	public:
		Particle2DBinnedOutput(int oi,int bcx,int bcy,int sc=0):interval(oi),
				binCountX(bcx),binCountY(bcy),stepCnt(sc),cartFout(0),gcFout(0) {
			cartParts.resize(binCountY);
			gcParts.resize(binCountY);
			
			if(binCountY<2) {
				cartFout=writeHeader("ParticleCartBinned.bin");
				gcFout=writeHeader("ParticleGCBinned.bin");

				fflush(cartFout);
				fflush(gcFout);
			}
		}
		
		~Particle2DBinnedOutput() {
			if(cartFout!=0) fclose(cartFout);
			if(gcFout!=0) fclose(gcFout);
		}

		template<class Population>
		void output(Population &pop) {
			if(stepCnt%interval==0) {
				for(int i=0; i<binCountY; ++i) {
					cartParts[i].clear();
					gcParts[i].clear();
				}
				float cartMiny=pop.gety(0);
				float cartMaxy=pop.gety(0);
				float gcMiny=pop.getY(0);
				float gcMaxy=pop.getY(0);
				for(int i=1; i<pop.getNumBodies(); ++i) {
					if(cartMiny>pop.gety(i)) cartMiny=pop.gety(i);
					if(cartMaxy<pop.gety(i)) cartMaxy=pop.gety(i);
					if(gcMiny>pop.getY(i)) gcMiny=pop.getY(i);
					if(gcMaxy<pop.getY(i)) gcMaxy=pop.getY(i);
				}
				for(int i=0; i<pop.getNumBodies(); ++i) {
					int cartBin=(int)((pop.gety(i)-cartMiny)*binCountY/(cartMaxy-cartMiny));
					int gcBin=(int)((pop.getY(i)-gcMiny)*binCountY/(gcMaxy-gcMiny));
					if(cartBin<0) cartBin=0;
						else if(cartBin>=binCountY) cartBin=binCountY-1;
					if(gcBin<0) gcBin=0;
						else if(gcBin>=binCountY) gcBin=binCountY-1;
					cartParts[cartBin].push_back(i);
					gcParts[gcBin].push_back(i);
				}
				CartCompareFunc<Population> cartComp(pop);
				for(int i=0; i<binCountY; ++i) {
					sort(cartParts[i].begin(),cartParts[i].end(),cartComp);
				}
				GCCompareFunc<Population> gcComp(pop);
				for(int i=0; i<binCountY; ++i) {
					sort(gcParts[i].begin(),gcParts[i].end(),gcComp);
				}
				if(binCountY>=2) {
					char buf[100];
					sprintf(buf,"ParticleCartBinned.%d.bin",stepCnt);
					cartFout=writeHeader(buf);
					sprintf(buf,"ParticleGCBinned.%d.bin",stepCnt);
					gcFout=writeHeader(buf);
				}
				for(int i=0; i<binCountY; ++i) {
					writeLineHeader(cartFout,cartMiny+(i+0.5)*(cartMaxy-cartMiny)/binCountY);
					writeLineHeader(gcFout,gcMiny+(i+0.5)*(gcMaxy-gcMiny)/binCountY);
					for(int j=0; j<binCountX; ++j) {
						PartBin bin;
						unsigned int start=j*cartParts[i].size()/binCountX;
						unsigned int end=(j+1)*cartParts[i].size()/binCountX;
						bin.minx=pop.getx(cartParts[i][start]);
						bin.maxx=pop.getx(cartParts[i][start]);
						for(unsigned int k=start; k<end && k<cartParts[i].size(); ++k) {
							bin.count++;
							double x=pop.getx(cartParts[i][k]);
							bin.x+=x;
							bin.y+=pop.gety(cartParts[i][k]);
							if(x<bin.minx) bin.minx=x;
							if(x>bin.maxx) bin.maxx=x;
							bin.tau+=3.14159*pop.getRadius(cartParts[i][k])*pop.getRadius(cartParts[i][k]);
							bin.e+=pop.gete(cartParts[i][k]);
							bin.i+=pop.geti(cartParts[i][k]);
							float phi=pop.getPhi(cartParts[i][k]);
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
							bin.tau*=binCountY/((bin.maxx-bin.minx)*(cartMaxy-cartMiny));
							bin.e/=bin.count;
							bin.i/=bin.count;
						}
						writeBin(cartFout,bin);
					}
					for(int j=0; j<binCountX; ++j) {
						PartBin bin;
						unsigned int start=j*gcParts[i].size()/binCountX;
						unsigned int end=(j+1)*gcParts[i].size()/binCountX;
						bin.minx=pop.getX(cartParts[i][start]);
						bin.maxx=pop.getX(cartParts[i][start]);
						for(unsigned int k=start; k<end && k<gcParts[i].size(); ++k) {
							bin.count++;
							double x=pop.getX(gcParts[i][k]);
							bin.x+=x;
							bin.y+=pop.getY(gcParts[i][k]);
							if(x<bin.minx) bin.minx=x;
							if(x>bin.maxx) bin.maxx=x;
							bin.tau+=3.14159*pop.getRadius(gcParts[i][k])*pop.getRadius(gcParts[i][k]);
							bin.e+=pop.gete(gcParts[i][k]);
							bin.i+=pop.geti(gcParts[i][k]);
							float phi=pop.getPhi(gcParts[i][k]);
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
							bin.tau*=binCountY/((bin.maxx-bin.minx)*(cartMaxy-cartMiny));
							bin.e/=bin.count;
							bin.i/=bin.count;
						}
						writeBin(gcFout,bin);
					}
				}
				if(binCountY>=2) {
					fclose(cartFout);
					fclose(gcFout);
					cartFout=0;
					gcFout=0;
				} else {
					fflush(cartFout);
					fflush(gcFout);
				}
			}
			stepCnt++;
		}
	private:
		FILE *writeHeader(const char *filename) {
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
		int stepCnt;
		
		vector<vector<int> > cartParts;
		vector<vector<int> > gcParts;

		FILE *cartFout;
		FILE *gcFout;
};

#endif
