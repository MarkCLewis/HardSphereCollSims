// BinaryDumpOutput.h
// This is an output routine that does a binary dump of the data to file.

#include <stdio.h>
#include <string>
//#include "GCPopulation.h"

template <class CoordType>
class BinaryDumpOutput {
	public:
		BinaryDumpOutput(int f):frequency(f),outCnt(0) {}

		BinaryDumpOutput(int f,int c):frequency(f),outCnt(c) {}

		template<class Population>
		void output(Population &pop) {
			if(outCnt%frequency==0) {
				char buf[100];
				sprintf(buf,"CartAndRad.%d.bin",outCnt);
				FILE *fout=fopen(buf,"wb");
				int num=pop.getNumBodies();
				fwrite(&num,sizeof(int),1,fout);
				fwrite(&(pop.getCart()[0]),sizeof(CoordType),num,fout);
				fwrite(&(pop.getRadius()[0]),sizeof(double),num,fout);
#ifdef SPIN
				fwrite(&(pop.getomega()[0]),sizeof(SpinVector),num,fout);
#endif
				fclose(fout);
			}
			outCnt++;
		}
	private:
		int frequency;
		int outCnt;
};
