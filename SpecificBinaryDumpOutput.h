// SpecificBinaryDumpOutput.h
// This is an output routine that does a binary dump of the data to file.

#include <stdio.h>
#include <string.h>
#include "GCPopulation.h"
#include "ProcessorCommunication.h"

#ifdef PARALLEL

class SpecificBinaryDumpOutput {
	public:
		SpecificBinaryDumpOutput(int f):frequency(f) {}

		bool doOnStep(int step) {
			return step%frequency==0;
		}

		void init(int step) {
			char buf[100];
			sprintf(buf,"CartAndRad.%d.bin",step);
			fout=fopen(buf,"wb");
			int num=0;
			radii.resize(0);

			// Write out a number just as a place holder.  We come
			// back and overwrite it.
			fwrite(&num,sizeof(int),1,fout);
		}

		void process(int num,vector<CartCoords> &c,vector<double> &r) {
			fwrite(&(c[0]),sizeof(CartCoords),num,fout);
			int oldSize=radii.size();
			radii.resize(oldSize+num);
			memcpy(&(radii[oldSize]),&(r[0]),num*sizeof(double));
		}

		void finalize(int step) {
			fwrite(&(radii[0]),sizeof(double),radii.size(),fout);

			fseek(fout,0,SEEK_SET);
			int sum=radii.size();
			fwrite(&sum,sizeof(int),1,fout);
			fclose(fout);
		}
	private:
		int frequency;
		std::vector<double> radii;
		FILE *fout;
};

#endif
