// TextOutput.h

#ifndef TEXT_OUTPUT
#define TEXT_OUTPUT

#include <stdio.h>

class TextOutput {
	public:
		TextOutput(int outputInterval):time(0.0),cnt(0),interval(outputInterval) {
			fout=fopen("ParticleText.txt","wt");
		}

		~TextOutput() {
			fclose(fout);
		}

		template<class Population>
		void output(Population &pop) {
			if(cnt%interval==0) {
			    for(int i=0; i<pop.getNumBodies(); i++) {
#ifndef SPIN
					fprintf(fout,"%e %e %e %e %e %e %e %e\n",cnt*pop.getTimeStep(),pop.getx(i),pop.gety(i),pop.getz(i),pop.getvx(i),pop.getvy(i),pop.getvz(i),pop.getRadius(i));
#else
					fprintf(fout,"%e %e %e %e %e %e %e %e %e %e %e\n",cnt*pop.getTimeStep(),pop.getx(i),pop.gety(i),pop.getz(i),pop.getvx(i),pop.getvy(i),pop.getvz(i),pop.getRadius(i),pop.getwx(i),pop.getwy(i),pop.getwz(i));
#endif
			    }
				fflush(fout);
			}
			cnt++;
		}
	private:
		FILE *fout;
		double time;
		int cnt;
		int interval;
};

#endif
	
