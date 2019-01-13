// DoubleForce.h
// This is the header file for the class that server as wrappers for
// the forcing.  These only need to be used if there are multiple forces
// that have to be applied.  That will probably be the case in all simulations
// however, because at the very least I need collisions and the satellite
// forcing.

template<class F1,class F2>
class DoubleForce {
	public:
		DoubleForce(F1 &f1,F2 &f2):force1(f1),force2(f2) {}

		template<class Population>
		void applyForce(Population &pop) {
			force1.applyForce(pop);
			force2.applyForce(pop);
		}
	private:
		F1 &force1;
		F2 &force2;
};
