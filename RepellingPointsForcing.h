/**
 * This class can be used to place a number of points that will repel particles.
 * The purpose of this is to help with the creation of non-uniform distributions.
 * Each point applies an outward force on nearby particles that has magnitude
 * mag at the point and linearly decays to zero at distance rad.
 * 
 * The first class is a forcing that can place a single repelling point.  The
 * second class places many and takes a vector of the first to initialize.
 */
 class RepellingPoint {
 	public:
 		RepellingPoint(double px,double py,double pz,double r,double m):
 			x(px),y(py),z(pz),rad(r),mag(m),cnt(0.0) {}
 			
 		template<class Population>
 		void applyForce(Population &pop) {
 			for(unsigned int i=0; i<pop.getNumBodies(); ++i) {
	 			double dx=pop.getx(i)-x;
	 			double dy=pop.gety(i)-y;
	 			double dz=pop.getz(i)-z;
	 			double dist=sqrt(dx*dx+dy*dy+dz*dz);
	 			if(dist<rad) {
	 				double tmag=(cnt<200.0)?(cnt/200.0):1.0;
	 				double fmag=tmag*pop.getTimeStep()*mag*(rad-dist)/rad/dist;
	 				pop.setvx(i,(1.0-0.2*fmag)*pop.getvx(i)+fmag*dx);
	 				pop.setvy(i,(1.0-0.2*fmag)*pop.getvy(i)+fmag*dy);
	 				pop.setvz(i,(1.0-0.2*fmag)*pop.getvz(i)+fmag*dz);
	 			}
 			}
 			cnt+=1.0;
 		}
 	private:
 		double x,y,z;
 		double rad;
 		double mag;
 		double cnt;
 };
 
class RepellingPointsForcing {
	public:
		RepellingPointsForcing(std::vector<RepellingPoint> &pnts):points(pnts) {}

		template<class Population>
		void applyForce(Population &pop) {
			for(unsigned int i=0; i<points.size(); ++i) {
				points[i].applyForce(pop);
			}
		}
	private:
		std::vector<RepellingPoint> points;
};
