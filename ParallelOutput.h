#ifndef PARALLEL_OUTPUT_H_
#define PARALLEL_OUTPUT_H_

#include <stdio.h>
#include <string.h>
#include "GCPopulation.h"
#include "ProcessorCommunication.h"

#ifdef PARALLEL

using std::vector;

template<class SpecificOutput>
class ParallelOutput {
	public:
		ParallelOutput(ProcessorCommunication &pc,SpecificOutput& so,int startCnt=0):procComm(pc),specificOutput(so),outCnt(startCnt) {
		}

		template<class Population>
		void output(Population &pop) {
			if(specificOutput.doOnStep(outCnt)) {
				if(procComm.getProcessNum()==0) {
					specificOutput.init(outCnt);
					specificOutput.process(pop.getNumReal(),pop.getCart(),pop.getRadius());

					// Loop through other processes.
					for(int i=1; i<procComm.getNumProcesses(); i++) {
						int tmp;
						numVect.resize(1);
						procComm.readFrom(i,numVect,tmp);
						int num=(int)numVect[0];
						buffer.resize(6*num);
						procComm.readFrom(i,buffer,tmp);
						radii.resize(num);
						procComm.readFrom(i,radii,tmp);
						
						// Process
						cart.resize(num);
						memcpy(&(cart[0]),&(buffer[0]),num*6*sizeof(double));
						specificOutput.process(num,cart,radii);
					}

					specificOutput.finalize(outCnt);
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
		ProcessorCommunication &procComm;
		SpecificOutput &specificOutput;
		int outCnt;

		vector<double> numVect;
		vector<double> buffer;
		vector<double> radii;
		vector<CartCoords> cart;
};

#endif

#endif 
