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
						ParticleIndex pi{i};
#ifndef SPIN
						fprintf(fout,"%e %e %e %e %e %e %e %e\n",cnt*pop.getTimeStep(),pop.getx(pi),pop.gety(pi),pop.getz(pi),pop.getvx(pi),pop.getvy(pi),pop.getvz(pi),pop.getRadius(pi));
#else
						fprintf(fout,"%e %e %e %e %e %e %e %e %e %e %e\n",cnt*pop.getTimeStep(),pop.getx(pi),pop.gety(pi),pop.getz(pi),pop.getvx(pi),pop.getvy(pi),pop.getvz(pi),pop.getRadius(pi),pop.getwx(pi),pop.getwy(pi),pop.getwz(pi));
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
	
