/*
 * CoefficientsOfRestitution.h
 *
 *  Created on: Jul 7, 2011
 *      Author: mlewis
 */

#ifndef COEFFICIENTSOFRESTITUTION_H_
#define COEFFICIENTSOFRESTITUTION_H_

class En_Bridges_Et_0p5 {
public:
	static double EpsilonN(double velCmPerSec) {
		double ret = 0.34*pow(velCmPerSec,-0.234);
		if(ret>1.0) return 1.0;
		else return ret;
	}
#ifdef SPIN
	static double EpsilonT(double velCmPerSec) {
		return 0.5;
	}
#endif
private:
};

template<int EtPercent>
class En_Bridges_Et {
public:
	static double EpsilonN(double velCmPerSec) {
		double ret = 0.34*pow(velCmPerSec,-0.234);
		if(ret>1.0) return 1.0;
		else return ret;
	}
#ifdef SPIN
	static double EpsilonT(double velCmPerSec) {
		return 0.01*EtPercent;
	}
#endif
private:
};

class En_0p5_Et_0p5 {
public:
	static double EpsilonN(double velCmPerSec) {
		if(velCmPerSec<1e-3) return 1.0;
		else return 0.5;
	}

#ifdef SPIN
	static double EpsilonT(double velCmPerSec) {
		return 0.5;
	}
#endif
private:
};

class En_0p25Bridges_Et_0p1 {
public:
	static double EpsilonN(double velCmPerSec) {
		double ret = 0.25*0.34*pow(velCmPerSec,-0.234);
		if(ret>1.0) return 1.0;
		else return ret;
	}
#ifdef SPIN
	static double EpsilonT(double velCmPerSec) {
		return 0.1;
	}
#endif
private:
};

#endif /* COEFFICIENTSOFRESTITUTION_H_ */
