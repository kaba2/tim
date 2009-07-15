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
	//% S is a real (m x n)-matrix that contains n samples of an
	//% m-dimensional signal.
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
	const mwSize samples = mxGetN(inputSet[0]);
	const mwSize dimension = mxGetM(inputSet[0]);

	real* rawData = mxGetPr(inputSet[0]);
	real maxRelativeError = *mxGetPr(inputSet[1]);
	integer kNearest = *mxGetPr(inputSet[2]);
	const integer threads = *mxGetPr(inputSet[3]);

	omp_set_num_threads(threads);

	const SignalPtr data = SignalPtr(
		new Signal(samples, dimension, rawData));
	
	outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	*rawResult = differentialEntropy(data, kNearest, 
		maxRelativeError, Euclidean_NormBijection<real>());
}
