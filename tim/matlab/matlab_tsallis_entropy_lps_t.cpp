// Description: tsallis_entropy_lps_t
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/tsallis_entropy_lps.h"

using namespace Tim;

namespace
{

	void matlabTemporalTsallisEntropyLps(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum
		{
			xIndex,
			timeWindowRadiusIndex,
			qIndex,
			kNearestSuggestionIndex,
			filterIndex,
			threadsIndex
		};

		std::vector<SignalPtr> xEnsemble;
		getSignals(inputSet[xIndex], std::back_inserter(xEnsemble));

		const integer timeWindowRadius = asInteger(inputSet[timeWindowRadiusIndex]);
		const real q = asReal(inputSet[qIndex]);
		const integer kNearestSuggestion = asInteger(inputSet[kNearestSuggestionIndex]);

		std::vector<real> filter;
		getReals(inputSet[filterIndex], std::back_inserter(filter));

		const integer threads = asInteger(inputSet[threadsIndex]);
		setNumberOfThreads(threads);

		const SignalPtr estimate = temporalTsallisEntropyLps(
			range(xEnsemble.begin(), xEnsemble.end()),
			timeWindowRadius, 
			q,
			kNearestSuggestion,
			range(filter.begin(), filter.end()));

		const integer nans = std::max(estimate->t(), 0);
		const integer skip = std::max(-estimate->t(), 0); 
		const integer samples = std::max(nans + estimate->samples() - skip, 0);

		outputSet[0] = mxCreateDoubleMatrix(1, samples, mxREAL);
		real* rawResult = mxGetPr(outputSet[0]);
		std::fill_n(rawResult, nans, nan<real>());
		std::copy(estimate->data().begin() + skip, 
			estimate->data().end(), rawResult + nans);
	}

	void addFunction()
	{
		matlabAddFunction(
			"tsallis_entropy_lps_t",
			matlabTemporalTsallisEntropyLps);
	}

	CallFunction run(addFunction);

}
