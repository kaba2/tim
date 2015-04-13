// Description: differential_entropy_sp
// DocumentationOf: differential_entropy_sp.m

#include "tim/matlab/tim_matlab.h"

#include "tim/core/differential_entropy_sp.h"

void force_linking_differential_entropy_sp() {};

using namespace Tim;

namespace
{

	void matlabDifferentialEntropySp(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			X,
			Inputs
		};

		enum Output
		{
			Estimate,
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		std::vector<Signal> xEnsemble = getSignals(inputSet[X]);

		real* outResult = matlabCreateScalar<real>(outputSet[Estimate]);
		*outResult = differentialEntropySp(
			countingRange(xEnsemble.begin(), xEnsemble.end()));
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_sp",
			matlabDifferentialEntropySp);
	}

	CallFunction run(addFunction);

}
