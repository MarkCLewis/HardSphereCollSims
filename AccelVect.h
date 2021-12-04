#ifndef ACCEL_VECT
#define ACCEL_VECT

struct AccelVect {
	AccelVect():ax(0.0),ay(0.0),az(0.0) {}
	AccelVect(double x, double y, double z): ax(x), ay(y), az(z) {}

	double ax,ay,az;

	void operator+=(const AccelVect &that) {
		ax += that.ax;
		ay += that.ay;
		az += that.az;
	}
};

#endif