// Coordinates.h
// This file defines the coordinate systems that are used.

#ifndef COORDINATES
#define COORDINATES

#include <math.h>

int step=0;

class CartCoords;

class GCCoords {
	public:
		static const double BETA;
		static const double A0;
		
		GCCoords() {}
		GCCoords(double vX,double vY,double ve,double vphi,double vi,double vzeta):
			X(vX),Y(vY),e(ve),phi(vphi),i(vi),zeta(vzeta) {}
		GCCoords(const CartCoords &cart) {
			set(cart);
		}
		GCCoords(const GCCoords &gc):X(gc.X),Y(gc.Y),e(gc.e),phi(gc.phi),
			i(gc.i),zeta(gc.zeta) { }

		double X,Y,e,phi,i,zeta; //,fill1,fill2;
		void set(const CartCoords &cart);
		void advance(double dt) {
			phi+=dt;
			zeta+=dt;
			Y-=2.0*A0*X*dt;
		}
		static double Ydot(double X) {
			return -2.0*A0*X;
		}
		static double phidot(double X) {
			return 1.0;
		}
		static double zetadot(double X) {
			return 1.0;
		}
	private:
};

const double GCCoords::BETA=2.0;
const double GCCoords::A0=0.75;

class EGCCoords {
	public:
		EGCCoords() {}
		EGCCoords(double vX,double vY,double ve,double vphi,double vi,double vzeta):
			X(vX),Y(vY),e(ve),phi(vphi),i(vi),zeta(vzeta) {}
		EGCCoords(const CartCoords &cart) {
			set(cart);
		}
		EGCCoords(const EGCCoords &gc):X(gc.X),Y(gc.Y),e(gc.e),phi(gc.phi),
			i(gc.i),zeta(gc.zeta) { }

		double X,Y,e,phi,i,zeta; //,fill1,fill2;
		void set(const CartCoords &cart);
		void advance(double dt) {
			double term=1.0-1.5*X;
			phi+=dt*term;
			zeta+=dt*term;
			Y-=dt*1.5*X*term;
		}
		static double Ydot(double X) {
			return -1.5*X*(1.0-1.5*X);
		}
		static double phidot(double X) {
			return 1.0-1.5*X;
		}
		static double zetadot(double X) {
			return 1.0-1.5*X;
		}
	private:
};

struct BasicCartCoords {
	double p[6];
};

class CartCoords {
	public:
		CartCoords() {}
		CartCoords(double px,double py,double pz,double pvx,double pvy,double pvz):
			x(px),y(py),z(pz),vx(pvx),vy(pvy),vz(pvz) {}
		CartCoords(const GCCoords &gc) {
			set(gc);
		}
		CartCoords(const EGCCoords &gc) {
			set(gc);
		}
		double x,y,z,vx,vy,vz; //,fill1,fill2;
		double get(int dim) {
			switch(dim) {
			case 0: return x;
			case 1: return y;
			case 2: return z;
			case 3: return vx;
			case 4: return vy;
			case 5: return vz;
//			case 6: return fill1;
//			case 7: return fill2;
			}
			printf("Bad dimension!!\n");
			return 0.0;
		}
		void set(const GCCoords &gc) {
			double cosp=cos(gc.phi),sinp=sin(gc.phi);
		
			x=gc.X-gc.e*cosp;
			y=gc.Y+GCCoords::BETA*gc.e*sinp;
			z=gc.i*cos(gc.zeta);
			vx=gc.e*sinp;
			vy=GCCoords::BETA*gc.e*cosp-2.0*GCCoords::A0*gc.X;
			vz=-gc.i*sin(gc.zeta);
		}
		void set(const EGCCoords &egc) {
			double cosp=cos(egc.phi),sinp=sin(egc.phi);
			x=egc.X-egc.e*cosp;
			y=egc.Y+GCCoords::BETA*egc.e*sinp;
			z=egc.i*cos(egc.zeta);
			vx=egc.e*sinp*(1-1.5*egc.X);
			vy=(1-1.5*egc.X)*(GCCoords::BETA*egc.e*cosp-2.0*GCCoords::A0*egc.X);
			vz=-egc.i*sin(egc.zeta);
		}
		double dist(const CartCoords &o) {
			double dx=x-o.x;
			double dy=y-o.y;
			double dz=z-o.z;
			return sqrt(dx*dx+dy*dy+dz*dz);
		}
	private:
};

void GCCoords::set(const CartCoords &cc) {
	X=(2.0*cc.x+cc.vy)*2.0;
	Y=cc.y-cc.vx*BETA;
	double dx=X-cc.x;
	double dy=cc.y-Y;
	e=sqrt(dx*dx+dy*dy/(BETA*BETA));
	phi=atan2(cc.y-Y,BETA*(X-cc.x));
	zeta=atan2(-cc.vz,cc.z);
	i=fabs(cc.z/cos(zeta));
}

void EGCCoords::set(const CartCoords &cc) {
	X=0.333333333*((1.0+6.0*cc.x)-sqrt((1.0+6.0*cc.x)*(1.0+6.0*cc.x)-6.0*(2.0*cc.vy+4.0*cc.x))); // !!! Consider doing a taylor expansion on the sqrt.
	double term=1.0-1.5*X;
	Y=cc.y-(cc.vx*GCCoords::BETA)/term;
	double esinp=cc.vx/term;
	double ecosp=(2.0*cc.vy+3.0*term*cc.x)/term;
	e=sqrt(esinp*esinp+ecosp*ecosp);
	phi=atan2(cc.vx,2.0*cc.vy+3.0*term*cc.x);
	zeta=atan2(-cc.vz,cc.z);
	i=fabs(cc.z/cos(zeta));
}

#ifdef SPIN
class SpinVector {
public:
	double x,y,z;
	SpinVector():x(0.0),y(0.0),z(0.0) {}
	SpinVector(double wx,double wy,double wz):x(wx),y(wy),z(wz) {}
};
#endif

#endif
