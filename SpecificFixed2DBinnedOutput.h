#ifndef SPECIFICFIXED2DBINNEDOUTPUT_H_
#define SPECIFICFIXED2DBINNEDOUTPUT_H_

#include <stdio.h>
#include <string.h>
#include "GCPopulation.h"
#include "ProcessorCommunication.h"

#ifdef PARALLEL

using std::vector;

struct PBin {
	PBin():count(0),gcCount(0),tau(0),e(0),i(0),vx(0),vy(0),gcTau(0),gcE(0),gcPhi(0) {}
	int count;
	int gcCount;
	float tau;
	float e;
	float i;
	float vx;
	float vy;
	float z;
	float gcTau;
	float gcE;
	float gcPhi;
	float gcI;
	float gcZeta;
};

template<class BoundaryCondition>
class SpecificFixed2DBinnedOutput {
	public:
		SpecificFixed2DBinnedOutput(int oi,int bcx,int bcy,double minX,double maxX,BoundaryCondition &bc):frequency(oi),binCountX(bcx),binCountY(bcy),minx(minX),maxx(maxX),bounds(bc) {
			printf("Do resizes to %d by %d\n",binCountY,binCountX);
			bins.resize(binCountY);
			for(int i=0; i<binCountY; ++i) {
				bins[i].resize(binCountX);
			}
			printf("Bins made\n");
			fflush(stdout);
		}

		bool doOnStep(int step) {
			return step%frequency==0;
		}

		void init(int step) {
			miny=bounds.getMinYTotal();
			maxy=bounds.getMaxYTotal();

			for(int i=0; i<binCountY; ++i) {
				for(int j=0; j<binCountX; ++j) {
					bins[i][j].count=0;
					bins[i][j].gcCount=0;
					bins[i][j].tau=0;
					bins[i][j].e=0;
					bins[i][j].i=0;
					bins[i][j].vx=0;
					bins[i][j].vy=0;
					bins[i][j].z=0;
					bins[i][j].gcTau=0;
					bins[i][j].gcE=0;
					bins[i][j].gcPhi=0;
					bins[i][j].gcI=0;
					bins[i][j].gcZeta=0;
				}
			}
		}

		void finalize(int step) {
			for(int i=0; i<binCountY; ++i) {
				for(int j=0; j<binCountX; ++j) {
					bins[i][j].tau*=(binCountX*binCountY)/((maxx-minx)*(maxy-miny));
					bins[i][j].gcTau*=(binCountX*binCountY)/((maxx-minx)*(maxy-miny));
					if(bins[i][j].count>0) {
						bins[i][j].e/=bins[i][j].count;
						bins[i][j].i/=bins[i][j].count;
						bins[i][j].vx/=bins[i][j].count;
						bins[i][j].vy/=bins[i][j].count;
						bins[i][j].z/=bins[i][j].count;
					}
					if(bins[i][j].gcCount>0) {
						bins[i][j].gcE/=bins[i][j].gcCount;
						bins[i][j].gcPhi/=bins[i][j].gcCount;
						bins[i][j].gcI/=bins[i][j].gcCount;
						bins[i][j].gcZeta/=bins[i][j].gcCount;
					}
				}
			}

			// Write to files.					
			char buf[100];
			sprintf(buf,"BinnedCart.%d.bin",step);
			writeToFile(buf,bins);
		}

		void process(int num,vector<CartCoords> &c,vector<double> &r) {
			for(int i=0; i<num; ++i) {
				GCCoords gc;
				gc.set(c[i]);
				while(gc.phi>3.14159) gc.phi-=2.0*3.14159;
				while(gc.phi<-3.14159) gc.phi+=2.0*3.14159;
				int binx=(int)((c[i].x-minx)*binCountX/(maxx-minx));
				int biny=(int)((c[i].y-miny)*binCountY/(maxy-miny));
				if(binx>=0 && binx<binCountX && biny>=0 && biny<binCountY) {
					bins[biny][binx].count++;
					bins[biny][binx].tau+=3.14159*r[i]*r[i];
					bins[biny][binx].e+=gc.e;
					bins[biny][binx].i+=gc.i;
					bins[biny][binx].vx+=c[i].vx;
					bins[biny][binx].vy+=c[i].vy;
					bins[biny][binx].z+=c[i].z;
				}
				binx=(int)((gc.X-minx)*binCountX/(maxx-minx));
				biny=(int)((gc.Y-miny)*binCountY/(maxy-miny));
				if(binx>=0 && binx<binCountX && biny>=0 && biny<binCountY) {
					bins[biny][binx].gcCount++;
					bins[biny][binx].gcTau+=3.14159*r[i]*r[i];
					bins[biny][binx].gcE+=gc.e;
					if(bins[biny][binx].count>1) {
						while(bins[biny][binx].gcPhi/(bins[biny][binx].gcCount-1)-gc.phi>3.14159) {
							gc.phi+=2.0*3.14159;
						} 
						while(bins[biny][binx].gcPhi/(bins[biny][binx].gcCount-1)-gc.phi<-3.14159) {
							gc.phi-=2.0*3.14159;
						}
					}
					bins[biny][binx].gcPhi+=gc.phi;
					bins[biny][binx].gcI+=gc.i;
					if(bins[biny][binx].count>1) {
						while(bins[biny][binx].gcZeta/(bins[biny][binx].gcCount-1)-gc.zeta>3.14159) {
							gc.zeta+=2.0*3.14159;
						}
						while(bins[biny][binx].gcZeta/(bins[biny][binx].gcCount-1)-gc.zeta<-3.14159) {
							gc.zeta-=2.0*3.14159;
						}
					}
					bins[biny][binx].gcZeta+=gc.zeta;
				}
			}
		}
		
	private:
		void writeToFile(char* filename,vector<vector<PBin> > &data) {
			FILE *fout=fopen(filename,"wb");

			// Write out header information.
			fwrite(&binCountX,sizeof(int),1,fout);
			fwrite(&binCountY,sizeof(int),1,fout);
			fwrite(&minx,sizeof(float),1,fout);
			fwrite(&maxx,sizeof(float),1,fout);
			fwrite(&miny,sizeof(float),1,fout);
			fwrite(&maxy,sizeof(float),1,fout);
			int numParams=2,numValues=11;
			fwrite(&numParams,sizeof(int),1,fout);
			fwrite("Count",sizeof(char),6,fout);
			fwrite("GC Count",sizeof(char),9,fout);
			fwrite(&numValues,sizeof(int),1,fout);
			fwrite("tau",sizeof(char),4,fout);
			fwrite("e",sizeof(char),2,fout);
			fwrite("i",sizeof(char),2,fout);
			fwrite("vx",sizeof(char),3,fout);
			fwrite("vy",sizeof(char),3,fout);
			fwrite("z",sizeof(char),2,fout);
			fwrite("GC tau",sizeof(char),7,fout);
			fwrite("GC e",sizeof(char),5,fout);
			fwrite("GC phi",sizeof(char),7,fout);
			fwrite("GC i",sizeof(char),5,fout);
			fwrite("GC zeta",sizeof(char),8,fout);

			for(unsigned int i=0; i<data.size(); ++i) {
				for(unsigned int j=0; j<data[i].size(); ++j) {
					fwrite(&(data[i][j].count),sizeof(int),1,fout);
					fwrite(&(data[i][j].gcCount),sizeof(int),1,fout);
					fwrite(&(data[i][j].tau),sizeof(float),1,fout);
					fwrite(&(data[i][j].e),sizeof(float),1,fout);
					fwrite(&(data[i][j].i),sizeof(float),1,fout);
					fwrite(&(data[i][j].vx),sizeof(float),1,fout);
					fwrite(&(data[i][j].vy),sizeof(float),1,fout);
					fwrite(&(data[i][j].z),sizeof(float),1,fout);
					fwrite(&(data[i][j].gcTau),sizeof(float),1,fout);
					fwrite(&(data[i][j].gcE),sizeof(float),1,fout);
					fwrite(&(data[i][j].gcPhi),sizeof(float),1,fout);
					fwrite(&(data[i][j].gcI),sizeof(float),1,fout);
					fwrite(&(data[i][j].gcZeta),sizeof(float),1,fout);
				}
			}
			fclose(fout);
		}
		
		int frequency;
		int binCountX;
		int binCountY;
		float minx,maxx,miny,maxy;
		vector<vector<PBin> > bins;
		BoundaryCondition &bounds;
};

#endif

#endif /*FIXED2DBINNEDOUTPUT_H_*/
