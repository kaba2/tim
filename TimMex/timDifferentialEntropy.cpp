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

	//% DIFFERENTIAL_ENTROPY
	//% A differential entropy estimate from samples.
	//%
	//% H = differential_entropy(S, epsilon, k, threads)
	//%
	//% where
	//%
	//% S is an arbitrary dimensional cell array whose linearization contains
	//% q trials of a signal. Each signal is a real (m x n)-matrix that 
	//% contains n samples of an m-dimensional signal.
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
	const mwSize trials = mxGetNumberOfElements(inputSet[0]);

	std::vector<SignalPtr> xEnsemble;
	xEnsemble.reserve(trials);

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

		xEnsemble.push_back(SignalPtr(
			new Signal(samples, dimension, rawData)));
	}

	const real maxRelativeError = *mxGetPr(inputSet[1]);
	const integer kNearest = *mxGetPr(inputSet[2]);
	const integer threads = *mxGetPr(inputSet[3]);

#if PASTEL_OMP != 0
	omp_set_num_threads(threads);
#endif
	
	outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	*rawResult = differentialEntropy(
		forwardRange(xEnsemble.begin(), xEnsemble.end()), 
		maxRelativeError,
		kNearest);
}
