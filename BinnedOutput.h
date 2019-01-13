// BinnedOutput.h
// This file has classes for binned output.  Right now I'm just writing a
// quick one to output optical depths so that I can look at my outputs and
// make sure that things are behaving well.  I think I will also want to add
// a class that does moving bins so that I can get the best resolution in
// areas of high particle counts.

#ifndef BINNED_OUTPUT
#define BINNED_OUTPUT

#include <stdio.h>

class BinnedOutput {
	public:
		BinnedOutput(int oi,int bd,int bc,int vd,double minv,double maxv):interval(oi),
				binDimension(bd),binCount(bc),valueDimension(vd),min(minv),max(maxv),stepCnt(0) {
			fout=fopen("density.txt","wt");
		}

		~BinnedOutput() {
			fclose(fout);
		}

		template<class Population>
		void output(Population &pop) {
			if(stepCnt%interval==0) {
				std::vector<double> sum(binCount);
				std::vector<double> cnt(binCount);
				for(int i=0; i<binCount; ++i) {
					sum[i]=0.0;
					cnt[i]=0.0;
				}
				for(int i=0; i<pop.getNumBodies(); ++i) {
					int bin=(int)((pop.get(i,binDimension)-min)*binCount/(max-min));
					if(bin>=0 && bin<binCount) {
						sum[bin]+=pop.get(i,valueDimension);
						cnt[bin]++;
					}
				}
				fprintf(fout,"%d ",stepCnt);
				for(int i=0; i<binCount; ++i) {
					if(cnt[i]==0) fprintf(fout,"0.0 ");
					else fprintf(fout,"%e ",sum[i]/cnt[i]);
				}
				fprintf(fout,"\n");
				fflush(fout);
			}
			stepCnt++;
		}
	private:
		FILE *fout;
		int interval;
		int binDimension;
		int binCount;
		int valueDimension;
		double min,max;
		int stepCnt;
};

#endif
