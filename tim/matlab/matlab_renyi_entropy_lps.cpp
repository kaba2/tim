// Description: renyi_entropy_lps
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/renyi_entropy_lps.h"

void force_linking_renyi_entropy_lps() {};

using namespace Tim;

namespace
{

	void matlabRenyiEntropyLps(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			X,
			Q,
			KNearestSuggestion,
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

		const real q = asReal(inputSet[Q]);
		const integer kNearestSuggestion = asInteger(inputSet[KNearestSuggestion]);

		outputSet[Estimate] = mxCreateDoubleMatrix(1, 1, mxREAL);
		real* rawResult = mxGetPr(outputSet[Estimate]);

		*rawResult = renyiEntropyLps(
			range(xEnsemble.begin(), xEnsemble.end()), 
			q, kNearestSuggestion);
	}

	void addFunction()
	{
		matlabAddFunction(
			"renyi_entropy_lps",
			matlabRenyiEntropyLps);
	}

	CallFunction run(addFunction);

}
