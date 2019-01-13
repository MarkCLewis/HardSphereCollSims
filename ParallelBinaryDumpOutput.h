// ParallelBinaryDumpOutput.h
// This is an output routine that does a binary dump of the data to file.

#include <stdio.h>
#include <string.h>
#include "GCPopulation.h"
#include "ProcessorCommunication.h"

#ifdef PARALLEL

class ParallelBinaryDumpOutput {
	public:
		ParallelBinaryDumpOutput(int f,ProcessorCommunication &pc):frequency(f),outCnt(0),procComm(pc) {}
		ParallelBinaryDumpOutput(int f,int startCnt,ProcessorCommunication &pc):frequency(f),outCnt(startCnt),procComm(pc) {}

		template<class Population>
		void output(Population &pop) {
			if(outCnt%frequency==0) {
				if(procComm.getProcessNum()==0) {
					char buf[100];
					sprintf(buf,"CartAndRad.%d.bin",outCnt);
					FILE *fout=fopen(buf,"wb");
					int num=pop.getNumReal();
					int sum=num;
					radii.resize(0);

					// Write out a number just as a place holder.  We come
					// back and overwrite it.
					fwrite(&num,sizeof(int),1,fout);
					fwrite(&(pop.getCart()[0]),sizeof(CartCoords),num,fout);

					// Loop through other processes.
					buffer.resize(1+7*pop.getNumReal());
					for(int i=1; i<procComm.getNumProcesses(); i++) {
						int tmp;
						buffer.resize(1);
						procComm.sendTo(i,buffer);
						procComm.readFrom(i,buffer,tmp);
						num=(int)buffer[0];
						sum+=num;
						printf("receive %d %d\n",num,sum);
						buffer.resize(6*num);
						procComm.readFrom(i,buffer,tmp);
						fwrite(&(buffer[0]),sizeof(CartCoords),num,fout);
						buffer.resize(num);
						procComm.readFrom(i,buffer,tmp);
						int oldSize=radii.size();
						radii.resize(radii.size()+num);
						memcpy(&(radii[oldSize]),&(buffer[0]),num*sizeof(double));
					}
					fwrite(&(pop.getRadius()[0]),sizeof(double),pop.getNumReal(),fout);
					fwrite(&(radii[0]),sizeof(double),radii.size(),fout);

					fseek(fout,0,SEEK_SET);
					fwrite(&sum,sizeof(int),1,fout);
					fclose(fout);
				} else {
					int tmp;
					buffer.resize(1);
					procComm.readFrom(0,buffer,tmp);
					buffer[0]=(double)pop.getNumReal();
					printf("About to send %d=%e\n",pop.getNumReal(),buffer[0]);
					procComm.sendTo(0,buffer);
					buffer.resize(6*pop.getNumReal());
					memcpy(&(buffer[0]),&(pop.getCart()[0]),pop.getNumReal()*6*sizeof(double));
					procComm.sendTo(0,buffer);
					buffer.resize(pop.getNumReal());
					memcpy(&(buffer[0]),&(pop.getRadius()[0]),pop.getNumReal()*sizeof(double));
					procComm.sendTo(0,buffer);
				}
			}
			outCnt++;
		}
	private:
		int frequency;
		int outCnt;
		ProcessorCommunication &procComm;
		std::vector<double> buffer;
		std::vector<double> radii;
};

#endif
