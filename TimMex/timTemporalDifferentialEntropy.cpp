#include "mex.h"

#include "tim/core/differential_entropy_kl.h"

#include <boost/static_assert.hpp>

#include <pastel/sys/pastelomp.h>

using namespace Tim;

void mexFunction(int outputs, mxArray *outputSet[],
				 int inputs, const mxArray *inputSet[])
{
	enum
	{
		RealIsDouble = boost::is_same<real, double>::value
	};
	BOOST_STATIC_ASSERT(RealIsDouble);

	//% TEMPORAL_DIFFERENTIAL_ENTROPY
	//% A temporal differential entropy estimate from samples.
	//%
	//% H = temporal_differential_entropy(S, timeWindowRadius, epsilon, k, threads)
	//%
	//% where
	//%
	//% S is an arbitrary dimensional cell array whose linearization contains
	//% q trials of a signal. Each signal is a real (m x n)-matrix that 
	//% contains n samples of an m-dimensional signal.
	//%
	//% TIMEWINDOWRADIUS determines the radius of the time-window in samples 
	//% inside which samples are taken into consideration to the estimate at 
	//% time instant t. This allows the estimate to be adaptive to temporal changes.
	//% If no such changes should happen, better accuracy can be 
	//% achieved by either setting 'timeWindowRadius' maximally wide
	//% or by using the differential_entropy() function instead.
	//%
	//% EPSILON is the maximum relative error in distance that
	//% nearest neighbor searching is allowed to result in.
	//% Higher tolerances result in enhanced performance, but
	//% increases errors in the estimate. Default 0.
	//%
	//% K determines which k:th nearest neighbor the algorithm
	//% uses for estimation. Default 1.
	//%
	//% THREADS determines the number of threads to use for parallelization.
	//% To fully take advantage of multiple cores in your machine, set this
	//% to the number of cores in your machine. Note however that this makes 
	//% your computer unresponsive to other tasks. When you need responsiveness, 
	//% spare one core for other work. Default 1 (no parallelization).

	// It is intentional to assign the width
	// and height the wrong way. The reason
	// is that Matlab uses column-major storage
	// while we use row-major storage.

	const mwSize signals = mxGetDimensions(inputSet[0])[0];
	const mwSize trials = mxGetDimensions(inputSet[0])[1];

	std::vector<SignalPtr> xEnsemble(trials);

	for (integer i = 0;i < trials;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[0], i);

		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.
		const mwSize samples = mxGetN(signalArray);
		const mwSize dimension = mxGetM(signalArray);

		real* rawData = mxGetPr(signalArray);

		xEnsemble[i] = SignalPtr(
			new Signal(samples, dimension, rawData));
	}

	const integer timeWindowRadius = *mxGetPr(inputSet[1]);
	const real maxRelativeError = *mxGetPr(inputSet[2]);
	const integer kNearest = *mxGetPr(inputSet[3]);
	const integer threads = *mxGetPr(inputSet[4]);

	omp_set_num_threads(threads);
	
	std::vector<real> estimateSet;

	temporalDifferentialEntropy(
		forwardRange(xEnsemble.begin(), xEnsemble.end()), 
		timeWindowRadius, std::back_inserter(estimateSet),
		kNearest, 
		maxRelativeError, Euclidean_NormBijection<real>());

	outputSet[0] = mxCreateDoubleMatrix(1, estimateSet.size(), mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	std::copy(estimateSet.begin(), estimateSet.end(),
		rawResult);
}
