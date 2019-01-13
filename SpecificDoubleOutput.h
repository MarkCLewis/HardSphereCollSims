// SpecificDoubleOutput.h
// This class is just like the DoubleForce class.  It is used simply so that
// I can bind together multiple output routines into a single class and have
// the templates hopefully inline the calls.  In theory, this class will never
// really appear in the code.

#ifndef SPECIFIC_DOUBLE_OUTPUT
#define SPECIFIC_DOUBLE_OUTPUT

#include<vector>

using std::vector;

template <class Output1,class Output2>
class SpecificDoubleOutput {
	public:
		SpecificDoubleOutput(Output1 &op1,Output2 &op2):o1(op1),o2(op2) { }

		bool doOnStep(int step) {
			do1=o1.doOnStep(step);
			do2=o2.doOnStep(step);
			return do1 || do2;
		}

		void init(int step) {
			if(do1) o1.init(step);
			if(do2) o2.init(step);
		}

		void process(int num,vector<CartCoords> &c,vector<double> &r) {
			if(do1) o1.process(num,c,r);
			if(do2) o2.process(num,c,r);
		}

		void finalize(int step) {
			if(do1) o1.finalize(step);
			if(do2) o2.finalize(step);
		}

	private:
		Output1 &o1;
		Output2 &o2;
		bool do1,do2;
};

#endif
