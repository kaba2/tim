#include "mex.h"

#include "tim/core/differential_entropy_kl.h"

#include <pastel/sys/pastelomp.h>

#include <boost/static_assert.hpp>

using namespace Tim;

void mexFunction(int outputs, mxArray *outputSet[],
				 int inputs, const mxArray *inputSet[])
{
	enum
	{
		RealIsDouble = boost::is_same<real, double>::value
	};
	BOOST_STATIC_ASSERT(RealIsDouble);

	enum
	{
		xIndex,
		timeWindowRadiusIndex,
		maxRelativeErrorIndex,
		kNearestIndex,
		threadsIndex
	};

	const integer trials = mxGetNumberOfElements(inputSet[xIndex]);

	std::vector<SignalPtr> xEnsemble;
	xEnsemble.reserve(trials);

	for (integer i = 0;i < trials;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[xIndex], i);

		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.
		const integer samples = mxGetN(signalArray);
		const integer dimension = mxGetM(signalArray);

		real* rawData = mxGetPr(signalArray);

		xEnsemble.push_back(SignalPtr(
			new Signal(samples, dimension, rawData)));
	}

	const integer timeWindowRadius = *mxGetPr(inputSet[timeWindowRadiusIndex]);
	const real maxRelativeError = *mxGetPr(inputSet[maxRelativeErrorIndex]);
	const integer kNearest = *mxGetPr(inputSet[kNearestIndex]);
	const integer threads = *mxGetPr(inputSet[threadsIndex]);

#if PASTEL_OMP != 0
	omp_set_num_threads(threads);
#endif
	
	std::vector<real> estimateSet;
	estimateSet.reserve(minSamples(forwardRange(xEnsemble.begin(), xEnsemble.end())));

	temporalDifferentialEntropy(
		forwardRange(xEnsemble.begin(), xEnsemble.end()), 
		timeWindowRadius, std::back_inserter(estimateSet),
		maxRelativeError,
		kNearest);

	outputSet[0] = mxCreateDoubleMatrix(1, estimateSet.size(), mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	std::copy(estimateSet.begin(), estimateSet.end(),
		rawResult);
}
