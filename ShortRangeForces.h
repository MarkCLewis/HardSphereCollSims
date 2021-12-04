/**
 * @file ShortRangeForces.cpp
 * @author mlewis@trinity.edu
 * @brief This file holds the classes that perform short-range forces.
 * @version 0.1
 * @date 2021-11-30
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef SHORT_RANGE_FORCES
#define SHORT_RANGE_FORCES

#include<vector>

#include "ParticleIndex.h"
#include "AccelVect.h"

using std::vector;

class ShortRangeGravityOnly {
  public:
    template<class Population>
    AccelVect applyForce(const Population &pop, ParticleIndex pi, ParticleIndex oi, double dx, double dy, double dz, double dist, double offsetX, double offsetY) {
      if (dist > pop.getRadius(pi) + pop.getRadius(oi)) {
        double mag = pop.getMass(oi) / (dist * dist * dist);
        return AccelVect(dx * mag, dy * mag, dz * mag);
      }
#ifndef SUPPRESS_OUT
      else if (dist < (pop.getRadius(pi) + pop.getRadius(oi))*0.99 && offsetX == 0.0 && offsetY == 0.0) {
        // Note that this can happen due to the initial conditions or the application of boundary conditions. 
        // It should not happen to more than a few particles per step after the first few steps.
        printf("Warning: Overlapping Gravity without offsets %d %d %e %e %e %f\n", pi.i, oi.i, dist, pop.getRadius(pi), pop.getRadius(oi), dist/(pop.getRadius(pi) + pop.getRadius(oi)));
      }
#endif
      return AccelVect(0.0, 0.0, 0.0);
    }
};

class ShortRangeGravitySpringDiscontinuous {
  public:
    ShortRangeGravitySpringDiscontinuous(double eps=0.5, double pdepth=0.02): epsilon(eps), penetrationDepth(pdepth) {
      logeps = log(epsilon);
      psqrlogesqr = logeps * logeps + 3.14159 * 3.14159;
      kconst = epsilon * (psqrlogesqr) / (3.14159 * 3.14159);
      printf("Constants %e %e %e\n", logeps, psqrlogesqr, kconst);
    }

    template<class Population>
    AccelVect applyForce(const Population &pop, ParticleIndex pi, ParticleIndex oi, double dx, double dy, double dz, double dist, double offsetX, double offsetY) {
      // The dx, dy, dz point from pi to oi.
      double mag = 0.0;
      double overlap = dist - (pop.getRadius(pi) + pop.getRadius(oi));
      if (overlap <= 0.0) {
        double mp = pop.getMass(pi);
        double mo = pop.getMass(oi);
        double mu = (mp * mo) / (mp + mo);
        double dr = penetrationDepth * (pop.getRadius(pi) + pop.getRadius(oi)) * 0.5;
        double v_i = (pop.getRadius(pi) + pop.getRadius(oi)) * 0.5;  // Assume that impact velocity scales as particle size.
        double k = mu * (v_i *v_i / (dr * dr)) * kconst;
        double c = 2 * logeps * sqrt(k * mu / psqrlogesqr);
        double vx = pop.getvx(pi) - pop.getvx(oi);
        double vy = pop.getvy(pi) - pop.getvy(oi);
        double vz = pop.getvz(pi) - pop.getvz(oi);
        double vxdx = vx * dx;
        double vydy = vy * dy;
        double vzdz = vz * dz;
        double vnorm = (vxdx + vydy + vzdz) / dist;
        mag = (k * std::max(overlap, -penetrationDepth * (pop.getRadius(pi) + pop.getRadius(oi))) + c * vnorm) / mp;
        return AccelVect(dx * mag / dist, dy * mag / dist, dz * mag / dist);
      } else {
        mag = pop.getMass(oi) / (dist * dist * dist);
        return AccelVect(dx * mag, dy * mag, dz * mag);
      }      
    }

  private:
    double epsilon;
    double penetrationDepth;
    double logeps;
    double psqrlogesqr;
    double kconst;
};

#endif // SHORT_RANGE_FORCES