#include "tim/matlab/tim_mex.h"

#include "tim/core/renyi_entropy_lps.h"

using namespace Tim;

void mexFunction(int outputs, mxArray *outputSet[],
				 int inputs, const mxArray *inputSet[])
{
	enum
	{
		xIndex,
		timeWindowRadiusIndex,
		qIndex,
		maxRelativeErrorIndex,
		kNearestIndex,
		threadsIndex
	};

	std::vector<SignalPtr> xEnsemble;
	getSignals(inputSet[xIndex], std::back_inserter(xEnsemble));

	const integer timeWindowRadius = getInteger(inputSet[timeWindowRadiusIndex]);
	const real q = getReal(inputSet[qIndex]);
	const real maxRelativeError = getReal(inputSet[maxRelativeErrorIndex]);
	const integer kNearest = getInteger(inputSet[kNearestIndex]);
	const integer threads = getInteger(inputSet[threadsIndex]);
	setNumberOfThreads(threads);

	std::vector<real> estimateSet;
	estimateSet.reserve(minSamples(forwardRange(xEnsemble.begin(), xEnsemble.end())));

	temporalRenyiEntropyLps(
		forwardRange(xEnsemble.begin(), xEnsemble.end()),
		timeWindowRadius, 
		std::back_inserter(estimateSet),
		q,
		maxRelativeError,
		kNearest);

	outputSet[0] = mxCreateDoubleMatrix(1, estimateSet.size(), mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	std::copy(estimateSet.begin(), estimateSet.end(),
		rawResult);
}
