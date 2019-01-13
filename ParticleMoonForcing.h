// ParticleMoonForcing.h
// This header file defines a type of forcing that can can be applied to a
// population of particles.  Instead of doing full particle to particle self-
// gravity, this class finds all particles above a certain size and does gravity
// between them and other particles.

class ParticleMoonForcing {
	public:
		ParticleMoonForcing(double rad):cutoffRadius(rad),count(0.0) {}
		
		template<class Population>
		void applyForce(Population &pop) {
			printf("Apply particle moon force.\n");
			fflush(stdout);
			bool forceDone=false;
			for(int i=0; i<pop.getNumBodies(); ++i) {
				if(pop.getRadius(i)>cutoffRadius) {
					forceDone=true;
					double massi=pop.getMass(i);
					if(count<100) massi*=count/100.0;
					for(int j=0; j<pop.getNumBodies(); ++j) {
						if(j!=i) {
							double dx=pop.getx(i)-pop.getx(j);
							double dy=pop.gety(i)-pop.gety(j);
							double dz=pop.getz(i)-pop.getz(j);
							double dist=sqrt(dx*dx+dy*dy+dz*dz);
							if(dist>pop.getRadius(i)+2.0*pop.getRadius(j)) {
								double massj=pop.getMass(j);
								double mag=pop.getTimeStep()*massj/(dist*dist*dist);
								pop.setvx(i,pop.getvx(i)-dx*mag);
								pop.setvy(i,pop.getvy(i)-dy*mag);
								pop.setvz(i,pop.getvz(i)-dz*mag);
								mag=pop.getTimeStep()*massi/(dist*dist*dist);
								pop.setvx(j,pop.getvx(j)+dx*mag);
								pop.setvy(j,pop.getvy(j)+dy*mag);
								pop.setvz(j,pop.getvz(j)+dz*mag);
							}
						}
					}
				}
			}
			count++;
			printf("Adjust particles.\n");
			fflush(stdout);
			if(forceDone) {
				for(int i=0; i<pop.getNumBodies(); ++i) {
					pop.adjustAfterForce(i);
				}
			}
		}
	private:
		double cutoffRadius;
		double count;
};

