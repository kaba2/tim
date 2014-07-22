// Description: divergence_wkv
// DocumentationOf: divergence_wkv.m

#include "tim/matlab/tim_matlab.h"

#include "tim/core/divergence_wkv.h"

void force_linking_divergence_wkv() {};

using namespace Tim;

namespace
{

	void matlabDivergenceWkv(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			X,
			Y,
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
		std::vector<Signal> yEnsemble = getSignals(inputSet[Y]);

		real* outResult = createScalar<real>(outputSet[Estimate]);
		*outResult = divergenceWkv(
			range(xEnsemble.begin(), xEnsemble.end()), 
			range(yEnsemble.begin(), yEnsemble.end()));
	}

	void addFunction()
	{
		matlabAddFunction(
			"divergence_wkv",
			matlabDivergenceWkv);
	}

	CallFunction run(addFunction);

}
