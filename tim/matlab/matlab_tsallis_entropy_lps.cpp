// Description: tsallis_entropy_lps
// DocumentationOf: tsallis_entropy_lps.m

#include "tim/matlab/tim_matlab.h"

#include "tim/core/tsallis_entropy_lps.h"

void force_linking_tsallis_entropy_lps() {};

using namespace Tim;

namespace
{

	void matlabTsallisEntropyLps(
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

		std::vector<Signal> xEnsemble = getSignals(inputSet[X]);

		real q = asScalar<real>(inputSet[Q]);
		integer kNearestSuggestion = asScalar<integer>(inputSet[KNearestSuggestion]);


		real* outResult = createScalar<real>(outputSet[Estimate]);
		*outResult = tsallisEntropyLps(
			countingRange(xEnsemble.begin(), xEnsemble.end()),
			q, kNearestSuggestion);
	}

	void addFunction()
	{
		matlabAddFunction(
			"tsallis_entropy_lps",
			matlabTsallisEntropyLps);
	}

	CallFunction run(addFunction);

}
