// DoubleOutput.h
// This class is just like the DoubleForce class.  It is used simply so that
// I can bind together multiple output routines into a single class and have
// the templates hopefully inline the calls.  In theory, this class will never
// really appear in the code.

template <class Output1,class Output2>
class DoubleOutput {
	public:
		DoubleOutput(Output1 &op1,Output2 &op2):o1(op1),o2(op2) { }

		template<class Population>
		void output(Population &pop) {
			o1.output(pop);
			o2.output(pop);
		}
	private:
		Output1 &o1;
		Output2 &o2;
};
