#ifndef FIXED_2D_BINNED_OUTPUT
#define FIXED_2D_BINNED_OUTPUT

#include <cstdio>
#include <vector>

using std::vector;
using std::fwrite;

struct Bin {
	Bin():count(0),gcCount(0),tau(0),e(0),i(0),vx(0),vy(0),gcTau(0),gcE(0),gcPhi(0) {}
	int count;
	int gcCount;
	float tau;
	float e;
	float i;
	float vx;
	float vy;
	float gcTau;
	float gcE;
	float gcPhi;
};

class Fixed2DBinnedOutput {
	public:
		Fixed2DBinnedOutput(int oi,int bcx,int bcy,double minX,double maxX,int sc=0):interval(oi),
				binCountX(bcx),binCountY(bcy),minx(minX),maxx(maxX),stepCnt(sc),fout(0) {
			bins.resize(binCountY);
			for(int i=0; i<binCountY; ++i) {
				bins[i].resize(binCountX);
			}
			
			if(binCountY<2) {
				writeHeader("FixedBinned.bin",0,0);
			}
		}
		
		~Fixed2DBinnedOutput() {
			if(fout!=0) fclose(fout);
		}

		template<class Population>
		void output(Population &pop) {
			if(stepCnt%interval==0) {
				for(int i=0; i<binCountY; ++i) {
					for(int j=0; j<binCountX; ++j) {
						bins[i][j].count=0;
						bins[i][j].gcCount=0;
						bins[i][j].tau=0;
						bins[i][j].e=0;
						bins[i][j].i=0;
						bins[i][j].vx=0;
						bins[i][j].vy=0;
						bins[i][j].gcTau=0;
						bins[i][j].gcE=0;
						bins[i][j].gcPhi=0;
					}
				}
				float miny=pop.gety(ParticleIndex{0});
				float maxy=pop.gety(ParticleIndex{0});
				for(ParticleIndex pi={0}; pi.i<pop.getNumBodies(); ++pi.i) {
					if(miny>pop.gety(pi)) miny=pop.gety(pi);
					if(maxy<pop.gety(pi)) maxy=pop.gety(pi);
				}
				for(ParticleIndex pi={0}; pi.i<pop.getNumBodies(); ++pi.i) {
					int binx=(int)((pop.getx(pi)-minx)*binCountX/(maxx-minx));
					int biny=(int)((pop.gety(pi)-miny)*binCountY/(maxy-miny));
					if(binCountY<2) biny=0;
					if(binx>=0 && binx<binCountX && biny>=0 && biny<binCountY) {
						bins[biny][binx].count++;
						bins[biny][binx].tau+=3.14159*pop.getRadius(pi)*pop.getRadius(pi);
						bins[biny][binx].e+=pop.gete(pi);
						bins[biny][binx].i+=pop.geti(pi);
						bins[biny][binx].vx+=pop.getvx(pi);
						bins[biny][binx].vy+=pop.getvy(pi);
					}
					binx=(int)((pop.getX(pi)-minx)*binCountX/(maxx-minx));
					biny=(int)((pop.getY(pi)-miny)*binCountY/(maxy-miny));
					if(binCountY<2) biny=0;
					if(binx>=0 && binx<binCountX && biny>=0 && biny<binCountY) {
						bins[biny][binx].gcCount++;
						bins[biny][binx].gcTau+=3.14159*pop.getRadius(pi)*pop.getRadius(pi);
						bins[biny][binx].gcE+=pop.gete(pi);
						float phi=pop.getPhi(pi);
						if(bins[biny][binx].count>1) {
							while(bins[biny][binx].gcPhi/(bins[biny][binx].gcCount-1)-phi>3.14159) {
								phi+=2.0*3.14159;
							} 
							while(bins[biny][binx].gcPhi/(bins[biny][binx].gcCount-1)-phi<-3.14159) {
								phi-=2.0*3.14159;
							}
						}
						bins[biny][binx].gcPhi+=phi;
					}
				}
				for(int i=0; i<binCountY; ++i) {
					for(int j=0; j<binCountX; ++j) {
						bins[i][j].tau*=(binCountX*binCountY)/((maxx-minx)*(maxy-miny));
						bins[i][j].gcTau*=(binCountX*binCountY)/((maxx-minx)*(maxy-miny));
						if(bins[i][j].count>0) {
							bins[i][j].e/=bins[i][j].count;
							bins[i][j].i/=bins[i][j].count;
							bins[i][j].vx/=bins[i][j].count;
							bins[i][j].vy/=bins[i][j].count;
						}
						if(bins[i][j].gcCount>0) {
							bins[i][j].gcE/=bins[i][j].gcCount;
							bins[i][j].gcPhi/=bins[i][j].gcCount;
						}
					}
				}
				
				if(binCountY>=2) {
					char buf[100];
					sprintf(buf,"FixedBinned.%d.bin",stepCnt);
					writeHeader(buf,miny,maxy);
					for(int i=0; i<binCountY; ++i) {
						writeLine(bins[i],0);
					}
					fclose(fout);
				} else {
					writeLine(bins[0],0.5*(miny+maxy));
				}
			}
			stepCnt++;
		}
	private:
		void writeHeader(const char* filename,float miny,float maxy) {
			fout=fopen(filename,"wb");

			// Write out header information.
			fwrite(&binCountX,sizeof(int),1,fout);
			fwrite(&binCountY,sizeof(int),1,fout);
			fwrite(&minx,sizeof(float),1,fout);
			fwrite(&maxx,sizeof(float),1,fout);
			fwrite(&miny,sizeof(float),1,fout);
			fwrite(&maxy,sizeof(float),1,fout);
			int numParams=2,numValues=8;
			fwrite(&numParams,sizeof(int),1,fout);
			fwrite("Count",sizeof(char),6,fout);
			fwrite("GC Count",sizeof(char),9,fout);
			fwrite(&numValues,sizeof(int),1,fout);
			fwrite("tau",sizeof(char),4,fout);
			fwrite("e",sizeof(char),2,fout);
			fwrite("i",sizeof(char),2,fout);
			fwrite("vx",sizeof(char),3,fout);
			fwrite("vy",sizeof(char),3,fout);
			fwrite("GC tau",sizeof(char),7,fout);
			fwrite("GC e",sizeof(char),5,fout);
			fwrite("GC phi",sizeof(char),7,fout);
			fflush(fout);
		}
	
		void writeLine(vector<Bin> &data,float y) {
			if(binCountY<2) fwrite(&y,sizeof(float),1,fout);
			for(unsigned int i=0; i<data.size(); ++i) {
				fwrite(&(data[i].count),sizeof(int),1,fout);
				fwrite(&(data[i].gcCount),sizeof(int),1,fout);
				fwrite(&(data[i].tau),sizeof(float),1,fout);
				fwrite(&(data[i].e),sizeof(float),1,fout);
				fwrite(&(data[i].i),sizeof(float),1,fout);
				fwrite(&(data[i].vx),sizeof(float),1,fout);
				fwrite(&(data[i].vy),sizeof(float),1,fout);
				fwrite(&(data[i].gcTau),sizeof(float),1,fout);
				fwrite(&(data[i].gcE),sizeof(float),1,fout);
				fwrite(&(data[i].gcPhi),sizeof(float),1,fout);
			}
			fflush(fout);
		}
		
		int interval;
		int binCountX;
		int binCountY;
		float minx,maxx;
		int stepCnt;
		
		vector<vector<Bin> > bins;
		
		FILE *fout;
};

#endif
