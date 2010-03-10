// Description: differential_entropy_kl_t
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/differential_entropy_kl.h"

using namespace Tim;

namespace
{

	void matlabTemporalDifferentialEntropyKl(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum
		{
			xIndex,
			timeWindowRadiusIndex,
			kNearestIndex,
			filterIndex,
			threadsIndex
		};

		std::vector<SignalPtr> xEnsemble;
		getSignals(inputSet[xIndex], std::back_inserter(xEnsemble));

		const integer timeWindowRadius = getInteger(inputSet[timeWindowRadiusIndex]);
		const integer kNearest = getInteger(inputSet[kNearestIndex]);

		std::vector<real> filter;
		getReals(inputSet[filterIndex], std::back_inserter(filter));

		const integer threads = getInteger(inputSet[threadsIndex]);
		setNumberOfThreads(threads);

		std::vector<real> estimateSet;
		estimateSet.reserve(minSamples(forwardRange(xEnsemble.begin(), xEnsemble.end())));

		temporalDifferentialEntropyKl(
			forwardRange(xEnsemble.begin(), xEnsemble.end()), 
			timeWindowRadius, std::back_inserter(estimateSet),
			kNearest,
			Default_NormBijection(),
			forwardRange(filter.begin(), filter.end()));

		outputSet[0] = mxCreateDoubleMatrix(1, estimateSet.size(), mxREAL);
		real* rawResult = mxGetPr(outputSet[0]);

		std::copy(estimateSet.begin(), estimateSet.end(),
			rawResult);
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_kl_t",
			matlabTemporalDifferentialEntropyKl);
	}

	CallFunction run(addFunction);

}
