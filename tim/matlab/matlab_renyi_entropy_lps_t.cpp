// Description: renyi_entropy_lps_t
// Documentation: tim_matlab_functions.txt

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

		const integer nans = std::max(estimate->t(), 0);
		const integer skip = std::max(-estimate->t(), 0); 
		const integer samples = std::max(nans + estimate->samples() - skip, 0);

		outputSet[Estimate] = mxCreateDoubleMatrix(1, samples, mxREAL);
		real* rawResult = mxGetPr(outputSet[Estimate]);
		std::fill_n(rawResult, nans, nan<real>());
		std::copy(estimate->data().begin() + skip, 
			estimate->data().end(), rawResult + nans);
	}

	void addFunction()
	{
		matlabAddFunction(
			"renyi_entropy_lps_t",
			matlabTemporalRenyiEntropyLps);
	}

	CallFunction run(addFunction);

}
