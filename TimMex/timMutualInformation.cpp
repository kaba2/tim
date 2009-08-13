#include "mex.h"

#include "tim/core/mutual_information_kraskow.h"

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

	//% MUTUAL_INFORMATION 
	//% A mutual information estimate from samples.
	//%
	//% I = mutual_information(X, Y, yLag, sigma, k, threads)
	//%
	//% where
	//%
	//% I is 1-dimensional row matrix which contains the mutual 
	//% information estimates [I(X, Y_1), ..., I(X, Y_p)],
	//% where Y_i is the signal y delayed by lagSet(i) samples.
	//%
	//% X is an arbitrary-dimensional cell-array whose linearization
	//% contains q trials of signal x.
	//%
	//% Y is an arbitrary-dimensional cell-array whose linearization 
	//% contains q trials of signal y.
	//%
	//% YLAG is the lag in samples which is applied to signal Y.
	//%
	//% SIGMA determines the radius of the time-window inside which
	//% samples are taken into consideration to the mutual information
	//% estimate at time instant t. The time window at time instant t
	//% is given by [t - sigma, t + sigma]. This allows the estimate to
	//% be adaptive to temporal changes in mutual information. If no
	//% such changes should happen, better accuracy can be achieved by
	//% setting 'sigma' larger than (half) the number of samples.
	//% Default: number of samples (i.e. everything).
	//%
	//% K determines which k:th nearest neighbor the algorithm
	//% uses for estimation. Default 1.
	//%
	//% THREADS determines the number of threads to use for parallelization.
	//% To fully take advantage of multiple cores in your machine, set this
	//% to the number of cores in your machine. Note however that this makes 
	//% your computer unresponsive to other tasks. When you need responsiveness, 
	//% spare one core for other work. Default 1 (no parallelization).
	//%
	//% Each signal is a real (m x n)-matrix that contains n samples of an
	//% m-dimensional signal. The dimension of X and Y need not coincide.
	//% However, the number of trials has to coincide.
	//% If the number of samples varies with trials, the function uses 
	//% the minimum sample count among the trials of X and Y.

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

	std::vector<SignalPtr> yEnsemble(trials);

	for (integer i = 0;i < trials;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[1], i);

		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.
		const mwSize samples = mxGetN(signalArray);
		const mwSize dimension = mxGetM(signalArray);

		real* rawData = mxGetPr(signalArray);

		yEnsemble[i] = SignalPtr(
			new Signal(samples, dimension, rawData));
	}

	const integer yLag = *mxGetPr(inputSet[2]);
	const integer sigma = *mxGetPr(inputSet[3]);
	const integer kNearest = *mxGetPr(inputSet[4]);
	const integer threads = *mxGetPr(inputSet[5]);

	omp_set_num_threads(threads);

	std::vector<real> result;
	mutualInformation(
		forwardRange(xEnsemble.begin(), xEnsemble.end()),
		forwardRange(yEnsemble.begin(), yEnsemble.end()),
		std::back_inserter(result),
		yLag,
		sigma,
		kNearest);

	outputSet[0] = mxCreateDoubleMatrix(1, result.size(), mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);
	std::copy(result.begin(), result.end(), rawResult);
}
