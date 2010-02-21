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
			maxRelativeErrorIndex,
			kNearestSuggestionIndex,
			threadsIndex
		};

		std::vector<SignalPtr> xEnsemble;
		getSignals(inputSet[xIndex], std::back_inserter(xEnsemble));

		const integer timeWindowRadius = getInteger(inputSet[timeWindowRadiusIndex]);
		const real q = getReal(inputSet[qIndex]);
		const real maxRelativeError = getReal(inputSet[maxRelativeErrorIndex]);
		const integer kNearestSuggestion = getInteger(inputSet[kNearestSuggestionIndex]);
		const integer threads = getInteger(inputSet[threadsIndex]);
		setNumberOfThreads(threads);

		std::vector<real> estimateSet;
		estimateSet.reserve(minSamples(forwardRange(xEnsemble.begin(), xEnsemble.end())));

		temporalTsallisEntropyLps(
			forwardRange(xEnsemble.begin(), xEnsemble.end()),
			timeWindowRadius, 
			std::back_inserter(estimateSet),
			q,
			maxRelativeError,
			kNearestSuggestion);

		outputSet[0] = mxCreateDoubleMatrix(1, estimateSet.size(), mxREAL);
		real* rawResult = mxGetPr(outputSet[0]);

		std::copy(estimateSet.begin(), estimateSet.end(),
			rawResult);
	}

	void addFunction()
	{
		matlabAddFunction(
			"tsallis_entropy_lps_t",
			matlabTemporalTsallisEntropyLps);
	}

	CallFunction run(addFunction);

}
