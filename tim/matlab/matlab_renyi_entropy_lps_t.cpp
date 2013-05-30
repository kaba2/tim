// Description: renyi_entropy_lps_t
// DocumentationOf: renyi_entropy_lps_t.m

#include "tim/matlab/tim_matlab.h"

#include "tim/core/renyi_entropy_lps.h"

void force_linking_renyi_entropy_lps_t() {};

using namespace Tim;

namespace
{

	void matlabTemporalRenyiEntropyLps(
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

		std::vector<SignalPtr> xEnsemble;
		getSignals(inputSet[X], std::back_inserter(xEnsemble));

		const integer timeWindowRadius = asScalar<integer>(inputSet[TimeWindowRadius]);
		const real q = asScalar<real>(inputSet[Q]);
		const integer kNearestSuggestion = asScalar<integer>(inputSet[KNearestSuggestion]);

		std::vector<real> filter;
		getScalars(inputSet[FilterIndex], std::back_inserter(filter));

		const SignalPtr estimate = temporalRenyiEntropyLps(
			range(xEnsemble.begin(), xEnsemble.end()),
			timeWindowRadius, 
			q,
			kNearestSuggestion,
			range(filter.begin(), filter.end()));

		const integer nans = std::max(estimate->t(), (integer)0);
		const integer skip = std::max(-estimate->t(), (integer)0); 
		const integer samples = std::max(nans + estimate->samples() - skip, (integer)0);

		Array<real> result = createArray<real>(samples, 1, 
			outputSet[Estimate]);
		std::fill_n(result.begin(), nans, nan<real>());
		std::copy(estimate->data().begin() + skip, 
			estimate->data().end(), result.begin() + nans);
	}

	void addFunction()
	{
		matlabAddFunction(
			"renyi_entropy_lps_t",
			matlabTemporalRenyiEntropyLps);
	}

	CallFunction run(addFunction);

}
