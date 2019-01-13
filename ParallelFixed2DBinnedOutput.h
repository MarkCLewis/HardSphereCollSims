#ifndef FIXED2DBINNEDOUTPUT_H_
#define FIXED2DBINNEDOUTPUT_H_

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
	float gcTau;
	float gcE;
	float gcPhi;
};

template<class BoundaryCondition>
class ParallelFixed2DBinnedOutput {
	public:
		ParallelFixed2DBinnedOutput(int oi,ProcessorCommunication &pc,int bcx,int bcy,double minX,double maxX,BoundaryCondition &bc,int startCnt=0):frequency(oi),
				binCountX(bcx),binCountY(bcy),minx(minX),maxx(maxX),outCnt(startCnt),procComm(pc),bounds(bc) {
			bins.resize(binCountY);
			for(int i=0; i<binCountY; ++i) {
				bins[i].resize(binCountX);
			}
		}

		template<class Population>
		void output(Population &pop) {
			if(outCnt%frequency==0) {
				float miny=bounds.getMinYTotal();
				float maxy=bounds.getMaxYTotal();
				if(procComm.getProcessNum()==0) {
					int num=pop.getNumReal();
					radii.resize(0);

					printf("cleared\n");
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

					// Process local particles.
					printf("Process local.\n");
					process(pop.getNumReal(),pop.getCart(),pop.getRadius(),miny,maxy);

					// Loop through other processes.
					for(int i=1; i<procComm.getNumProcesses(); i++) {
						printf("Process %d\n",i);
						int tmp;
						numVect.resize(1);
						procComm.readFrom(i,numVect,tmp);
						num=(int)numVect[0];
						printf("Resize buffer %d\n",num);
						buffer.resize(6*num);
						printf("Reading buffer.\n");
						procComm.readFrom(i,buffer,tmp);
						printf("Resizing radii.\n");
						radii.resize(num);
						printf("Reading radii.\n");
						procComm.readFrom(i,radii,tmp);
						
						// Process
						printf("Resize cart.\n");
						cart.resize(num);
						printf("Memcpy\n");
						memcpy(&(cart[0]),&(buffer[0]),pop.getNumReal()*6*sizeof(double));
						printf("Do process\n");
						process(num,cart,radii,miny,maxy);
					}

					printf("Finalize bins\n");
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

					// Write to files.					
					printf("Output file\n");
					char buf[100];
					sprintf(buf,"BinnedCart.%d.bin",outCnt);
					writeToFile(buf,bins,miny,maxy);
					printf("Done with master binning.\n");
				} else {
					numVect.resize(1);
					numVect[0]=(double)pop.getNumReal();
					procComm.sendTo(0,numVect);
					buffer.resize(6*pop.getNumReal());
					memcpy(&(buffer[0]),&(pop.getCart()[0]),pop.getNumReal()*6*sizeof(double));
					procComm.sendTo(0,buffer);
					radii.resize(pop.getNumReal());
					memcpy(&(radii[0]),&(pop.getRadius()[0]),pop.getNumReal()*sizeof(double));
					procComm.sendTo(0,radii);
				}
			}
			outCnt++;
		}
	private:
		void process(int num,vector<CartCoords> &c,vector<double> &r,float miny,float maxy) {
			for(unsigned int i=0; i<c.size(); ++i) {
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
				}
			}
		}
		
		void writeToFile(char* filename,vector<vector<PBin> > &data,float miny,float maxy) {
			FILE *fout=fopen(filename,"wb");

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

			for(unsigned int i=0; i<data.size(); ++i) {
				for(unsigned int j=0; j<data[i].size(); ++j) {
					fwrite(&(data[i][j].count),sizeof(int),1,fout);
					fwrite(&(data[i][j].gcCount),sizeof(int),1,fout);
					fwrite(&(data[i][j].tau),sizeof(float),1,fout);
					fwrite(&(data[i][j].e),sizeof(float),1,fout);
					fwrite(&(data[i][j].i),sizeof(float),1,fout);
					fwrite(&(data[i][j].vx),sizeof(float),1,fout);
					fwrite(&(data[i][j].vy),sizeof(float),1,fout);
					fwrite(&(data[i][j].gcTau),sizeof(float),1,fout);
					fwrite(&(data[i][j].gcE),sizeof(float),1,fout);
					fwrite(&(data[i][j].gcPhi),sizeof(float),1,fout);
				}
			}
			fclose(fout);
		}
		
		int frequency;
		int binCountX;
		int binCountY;
		float minx,maxx;
		int outCnt;
		ProcessorCommunication &procComm;
		vector<double> numVect;
		vector<double> buffer;
		vector<double> radii;
		vector<CartCoords> cart;
		vector<vector<PBin> > bins;
		BoundaryCondition &bounds;
};

#endif

#endif /*FIXED2DBINNEDOUTPUT_H_*/
