// Description: tsallis_entropy_lps_t
// DocumentationOf: tsallis_entropy_lps_t.m

#include "tim/corematlab/tim_matlab.h"

#include "tim/core/tsallis_entropy_lps.h"

void force_linking_tsallis_entropy_lps_t() {};

using namespace Tim;

namespace
{

	void matlabTemporalTsallisEntropyLps(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			X,
			TimeWindowRadius,
			Q,
			KNearestSuggestion,
			FilterIndex,
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

		integer timeWindowRadius = matlabAsScalar<integer>(inputSet[TimeWindowRadius]);
		real q = matlabAsScalar<real>(inputSet[Q]);
		integer kNearestSuggestion = matlabAsScalar<integer>(inputSet[KNearestSuggestion]);

		std::vector<real> filter;
		matlabGetScalars(inputSet[FilterIndex], std::back_inserter(filter));

		Signal estimate = temporalTsallisEntropyLps(
			countingRange(xEnsemble.begin(), xEnsemble.end()),
			timeWindowRadius, 
			q,
			kNearestSuggestion,
			range(filter.begin(), filter.end()));

		integer nans = std::max(estimate.t(), (integer)0);
		integer skip = std::max(-estimate.t(), (integer)0); 
		integer samples = std::max(nans + estimate.samples() - skip, (integer)0);

		Array<real> result = matlabCreateArray<real>(
			samples, 1, outputSet[Estimate]);
		std::fill_n(result.begin(), nans, (real)Nan());
		std::copy(estimate.data().begin() + skip, 
			estimate.data().end(), result.begin() + nans);
	}

	void addFunction()
	{
		matlabAddFunction(
			"tsallis_entropy_lps_t",
			matlabTemporalTsallisEntropyLps);
	}

	CallFunction run(addFunction);

}
