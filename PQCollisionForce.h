#ifndef PQCollisionForce
#define PQCollisionForce

#include<vector>
#include<math.h>

#include "ParticleIndex.h"
#include "AccelVect.h"

using std::vector;

struct PQForceEvent {
  double time;
  ParticleIndex p1;
  ParticleIndex p2;
};

class ShortRangeGravitySpringDiscontinuousWithPQ {
  public:
    ShortRangeGravitySpringDiscontinuousWithPQ(vector<PQForceEvent> &vect, double eps=0.5, double pdepth=0.02): 
      pq(vect),
      epsilon(eps), 
      penetrationDepth(pdepth), 
      logeps(log(epsilon)), 
      psqrlogesqr(logeps * logeps + 3.14159 * 3.14159), 
      kconst(epsilon * (psqrlogesqr) / (3.14159 * 3.14159)) {}

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
#ifndef SUPPRESS_OUT        
        if (overlap < -20 * dr) {
          printf("Warning: large soft-sphere overlap: %d %d %3.2f%% %e %e\n", pi.i, oi.i, -overlap / (pop.getRadius(pi) + pop.getRadius(oi)) * 100.0, pop.getRadius(pi), pop.getRadius(oi));
        }
#endif
        mag = (k * overlap + c * vnorm) / mp;

        return AccelVect(dx * mag / dist, dy * mag / dist, dz * mag / dist);
      } else {
        mag = pop.getMass(oi) / (dist * dist * dist);
        return AccelVect(dx * mag, dy * mag, dz * mag);
      }      
    }

  private:
    const double epsilon;
    const double penetrationDepth;
    const double logeps;
    const double psqrlogesqr;
    const double kconst;
    const vector<PQForceEvent> &pq;
};

class PQCollisionForce {

};

#endif