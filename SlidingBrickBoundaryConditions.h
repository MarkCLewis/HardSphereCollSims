// BoundaryConditions.h
// This file contains some classes for different boundary conditions.  These
// classes are one of the main places where code for multiple processors shows
// up, because without gravity using boundaing layers the collisions run like
// they are on a single machine for each timestep.  At the end of the timestep
// the talk to one another about which particles have moved across the lines.
// A significant question here is whether I should make separate versions for
// the sinlge processor BCs or if that should be special code in a more
// complete version.

/*
Modification Log:
  July 28, 2006 - added implementation for sliding brick boundary conditions with 8 mirrors by ...
                  - added shearOffset variable to the class's private stock
                  - added shearOffset(0.) to class constructor
                  - added shearOffset calculation that iteratively subtracts 2*A0 times the cell height times the step size
                  - added method getShearOffset() to return the shear offset when called by GravTreeMirrors.h

*/
#ifndef BOUNDARY_CONDITIONS
#define BOUNDARY_CONDITIONS

#include <vector>
#include "GCPopulation.h"
#include "ProcessorCommunication.h"
#include "Coordinates.h"

class SlidingBrick {
	public:
		SlidingBrick(double minX,double maxX,double minY,double maxY):minx(minX),maxx(maxX),miny(minY),maxy(maxY),shearOffset(0.) {}

		template<class Population>
		void apply(Population &pop) {
			#pragma omp parallel for
			for(int ii = 0; ii<pop.getNumBodies(); ++ii) {
				ParticleIndex pi={ii};
				while(pop.getx(pi)<minx) {
					pop.setx(pi,pop.getx(pi)+(maxx-minx));
					double newy=pop.gety(pi)+2*shearOffset;
					if(newy>maxy) newy-=(maxy-miny);
					if(newy<miny) newy+=(maxy-miny);
					pop.sety(pi,newy);
					pop.setvy(pi,pop.getvy(pi)-2.0*GCCoords::A0*(maxx-minx));
					pop.adjustAfterForce(pi);
				} 
				while(pop.getx(pi)>maxx) {
					pop.setx(pi,pop.getx(pi)-(maxx-minx));
					double newy=pop.gety(pi)-2*shearOffset;
					if(newy>maxy) newy-=(maxy-miny);
					if(newy<miny) newy+=(maxy-miny);
					pop.sety(pi,newy);
					pop.setvy(pi,pop.getvy(pi)+2.0*GCCoords::A0*(maxx-minx));
					pop.adjustAfterForce(pi);
				}
				while(pop.gety(pi)<miny) {
					pop.sety(pi,pop.gety(pi)+(maxy-miny));
					pop.adjustAfterForce(pi);
				} 
				while(pop.gety(pi)>maxy) {
					pop.sety(pi,pop.gety(pi)-(maxy-miny));
					pop.adjustAfterForce(pi);
				}
			}
            
      shearOffset -= 2*GCCoords::A0*(maxx-minx)*0.5*pop.getTimeStep();
      if(shearOffset < miny) { shearOffset += maxy-miny; }
		}

		bool checkIfIn(double x,double y) {
			return (x>=minx && x<=maxx && y>=miny && y<=maxy);
		}

		double getArea() {
			return (maxx-minx)*(maxy-miny);
		}

		double getMinX() { return minx; }
		double getMaxX() { return maxx; }
		double getMinY() { return miny; }
		double getMaxY() { return maxy; }
    
		double getShearOffset() { return shearOffset;}

	private:
		double minx,maxx,miny,maxy;
		double shearOffset;
};

#endif
