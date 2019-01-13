/*
 * ConstantAccel.h
 *
 *  Created on: Jun 19, 2010
 *      Author: mlewis
 */

#ifndef CONSTANTACCEL_H_
#define CONSTANTACCEL_H_

class ConstantAccel {
	public:
		ConstantAccel(double a=9.8):g(a) {}

		template <class Population>
		void applyForce(Population &pop) {
			int nb=pop.getNumBodies();
			double step=pop.getTimeStep();
			for(int i=0; i<nb; ++i) {
				pop.setvz(i,pop.getvz(i)-g*step);
			}
		}
	private:
		double g;
};

#endif /* CONSTANTACCEL_H_ */
