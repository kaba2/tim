#include "mex.h"

#include "tim/core/partial_mutual_information.h"

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

	//% PARTIAL_MUTUAL_INFORMATION 
	//% A partial mutual information estimate from samples.
	//%
	//% I = partial_mutual_information(X, Y, Z, yLag, zLag, k, threads)
	//%
	//% where
	//%
	//% X, Y, and Z are arbitrary-dimensional cell-arrays whose 
	//% linearizations contain q trials of signal x, y, and z, 
	//% respectively.
	//%
	//% YLAG and ZLAG are the lags in samples applied to signals
	//% y and z, respectively.
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
	//% m-dimensional signal. The dimensions of X, Y, and Z need not coincide.
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

	std::vector<SignalPtr> zEnsemble(trials);

	for (integer i = 0;i < trials;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[2], i);

		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.
		const mwSize samples = mxGetN(signalArray);
		const mwSize dimension = mxGetM(signalArray);

		real* rawData = mxGetPr(signalArray);

		zEnsemble[i] = SignalPtr(
			new Signal(samples, dimension, rawData));
	}

	const integer yLag = *mxGetPr(inputSet[3]);
	const integer zLag = *mxGetPr(inputSet[4]);
	const integer kNearest = *mxGetPr(inputSet[5]);
	const integer threads = *mxGetPr(inputSet[6]);

	omp_set_num_threads(threads);

	outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	*rawResult = partialMutualInformation(
		forwardRange(xEnsemble.begin(), xEnsemble.end()),
		forwardRange(yEnsemble.begin(), yEnsemble.end()),
		forwardRange(zEnsemble.begin(), zEnsemble.end()),
		yLag,
		zLag,
		kNearest);
}
