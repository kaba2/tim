// Description: differential_entropy_kl
// DocumentationOf: differential_entropy_kl.m

#include "tim/matlab/tim_matlab.h"

#include "tim/core/differential_entropy_kl.h"

void force_linking_differential_entropy_kl() {};

using namespace Tim;

namespace
{

	void matlabDifferentialEntropyKl(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			X,
			KNearest,
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

		integer kNearest = matlabAsScalar<integer>(inputSet[KNearest]);


		real* outResult = matlabCreateScalar<real>(outputSet[Estimate]);
		*outResult = differentialEntropyKl(
			countingRange(xEnsemble.begin(), xEnsemble.end()), 
			kNearest);
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_kl",
			matlabDifferentialEntropyKl);
	}

	CallFunction run(addFunction);

}
