// Description: divergence_wkv
// Documentation: tim_matlab_functions.txt

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

		std::vector<SignalPtr> xEnsemble;
		getSignals(inputSet[X], std::back_inserter(xEnsemble));

		std::vector<SignalPtr> yEnsemble;
		getSignals(inputSet[Y], std::back_inserter(yEnsemble));

		outputSet[Estimate] = mxCreateDoubleMatrix(1, 1, mxREAL);
		real* rawResult = mxGetPr(outputSet[Estimate]);

		*rawResult = divergenceWkv(
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
