#include<cstdio>
#include<vector>
#include<cmath>

#include "Coordinates.h"

struct vect3 {
  double x, y, z;
  vect3(): x(0.0), y(0.0), z(0.0) {}
  vect3(double x, double y, double z): x(x), y(y), z(z) {}
  vect3 operator+(const vect3 &that) const { return vect3(x+that.x, y+that.y, z+that.z); }
  vect3 operator-(const vect3 &that) const { return vect3(x-that.x, y-that.y, z-that.z); }
  vect3 operator*(double c) const { return vect3(x*c, y*c, z*c); }
  vect3 operator/(double c) const { return vect3(x/c, y/c, z/c); }
  vect3 operator+=(const vect3 &that) { 
    x += that.x; y += that.y; z += that.z;
    return *this;
  }
  vect3 operator-=(const vect3 &that) { 
    x -= that.x; y -= that.y; z -= that.z;
    return *this;
  }
  vect3 operator*=(double c) { 
    x += c; y *= c; z *= c;
    return *this;
  }
  vect3 operator/=(double c) { 
    x /= c; y /= c; z /= c;
    return *this;
  }
  double dot(const vect3 &that) const { return x*that.x + y*that.y + z*that.z; }
};

using std::vector;

inline double sigmoid(double x) { 
  return 1.0/(1.0+exp(-x)); 
}
inline double sigmoidP(double x) {
  if (x > 100.0) return 0.0;
  double denom = exp(x)+1.0;
  return exp(x)/(denom * denom); 
}

void calcAccJerk(const vector<vect3> &pos, const vector<vect3> &vel, const vector<double> &rad, double rho, 
        vector<vect3> &acc, vector<vect3> &jerk) {
  for (vector<CartCoords>::size_type i = 0; i < pos.size(); ++i) {
    acc[i] = {0.0, 0.0, 0.0};
    jerk[i] = {0.0, 0.0, 0.0};
  }
  for (vector<CartCoords>::size_type i = 0; i < pos.size(); ++i) {
    for (vector<CartCoords>::size_type j = i+1; j < pos.size(); ++j) {
      vect3 rji = pos[j] - pos[i];
      vect3 vji = vel[j] - vel[i];
      double r2 = rji.dot(rji);
      double rv_r2 = rji.dot(vji) / r2;
      double r = sqrt(r2);
      double r3 = r * r2;

      vect3 da = rji / r3;
      vect3 dj = (vji - rji * 3 * rv_r2) / r3;
      double massi = 4 * M_PI * rho / 3 * rad[i] * rad[i] * rad[i];
      double massj = 4 * M_PI * rho / 3 * rad[j] * rad[j] * rad[j];
      acc[i] += da * massj;
      acc[j] -= da * massi;
      jerk[i] += dj * massj;
      jerk[j] -= dj * massi;
    }
  }
}

void predictStep(vector<vect3> &pos, vector<vect3> &vel,
        const vector<vect3> &acc, const vector<vect3> &jerk, double dt) {
  double dt2 = dt*dt/2.0;
  double dt3 = dt2*dt/6.0;

  for (vector<CartCoords>::size_type i = 0; i < pos.size(); ++i) {
    // TODO: These might not be optimally efficient without expression templates.
    pos[i] += vel[i] * dt + acc[i] * dt2 + jerk[i] * dt3;
    vel[i] += acc[i] * dt + jerk[i] * dt2;
  }
}

void correctStep(vector<vect3> &pos, vector<vect3> &vel, const vector<vect3> &acc, const vector<vect3> &jerk, 
        const vector<vect3> &oldPos, const vector<vect3> &oldVel, vector<vect3> &oldAcc, vector<vect3> &oldJerk, double dt) {
  // TODO: These might not be optimally efficient without expression templates.
  for (vector<CartCoords>::size_type i = 0; i < pos.size(); ++i) {
    vel[i] = oldVel[i] + (oldAcc[i] + acc[i])*dt/2 + (oldJerk[i] - jerk[i])*dt*dt/12;
    pos[i] = oldPos[i] + (oldVel[i] + vel[i])*dt/2 + (oldAcc[i] - acc[i])*dt*dt/12;
  }
}

void evolveStep(vector<vect3> &pos,vector<vect3> &vel, const vector<double> &rad, const double rho, 
        vector<vect3> &acc, vector<vect3> &jerk, double dt) {
  // TODO: This isn't ideally efficient. Better to keep two vectors around and reuse, but it will do for now.
  vector<vect3> oldPos = {pos};
  vector<vect3> oldVel = {vel};
  vector<vect3> oldAcc = {acc};
  vector<vect3> oldJerk = {jerk};

  predictStep(pos, vel, acc, jerk, dt);
  calcAccJerk(pos, vel, rad, rho, acc, jerk);
  correctStep(pos, vel, acc, jerk, oldPos, oldVel, oldAcc, oldJerk, dt);
}

const double n = 1.0;
const double kappa = 1.0;
const double n_z = 1.0;

// Call this after all others so the accels are completed by the time this is called.
void hillsForce(const vect3 &pos, const vect3 &vel, vect3 &acc, vect3 &jerk) {
  acc.x += 2.0 * n * vel.y - (kappa * kappa - 4.0 * n * n) * pos.x;
  acc.y += -2.0 * n * vel.x;
  acc.z += -n_z * n_z * pos.z;
  jerk.x += 2.0 * n * acc.y - (kappa * kappa - 4.0 * n * n) * vel.x;
  jerk.y += -2.0 * n * acc.x;
  jerk.z += -n_z * n_z * vel.z;
}

int main(int argc, char **argv) {
  double dt = 0.001 * 2 *M_PI;
  vector<vect3> pos = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
  vector<vect3> vel = {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
  const vector<double> rad = {1e-7, 1e-14};
  const double rho = 3.0 / (4.0 * M_PI *rad[0]*rad[0]*rad[0]);
  vector<vect3> acc(pos.size());
  vector<vect3> jerk(pos.size());

  calcAccJerk(pos, vel, rad, rho, acc, jerk);
  for (double t = 0.0; t < 2e5 * M_PI; t += dt) {
    evolveStep(pos, vel, rad, rho, acc, jerk, dt);
  }
  printf("dist = %e\n", sqrt((pos[0]-pos[1]).dot(pos[0]-pos[1])));
  return 0;
}