#include "mex.h"

#include "tim/core/transfer_entropy.h"

#include <pastel/sys/pastelomp.h>

#include <boost/static_assert.hpp>

#include <algorithm>

using namespace Tim;

void mexFunction(int outputs, mxArray *outputSet[],
				 int inputs, const mxArray *inputSet[])
{
	enum
	{
		RealIsDouble = boost::is_same<real, double>::value
	};
	BOOST_STATIC_ASSERT(RealIsDouble);

	//% TRANSFER_ENTROPY 
	//% A multivariate transfer entropy estimate from samples.
	//%
	//% I = transfer_entropy(X, W, Y, Z, sigma, k, threads)
	//%
	//% where
	//%
	//% X is a cell-array of arbitrary dimension which contains q trials of 
	//% the already embedded X-signal. Call the signal from which the X-signal
	//% was embedded from the x-signal.
	//%
	//% W is a cell-array of arbitrary dimension which contains q trials of
	//% the time-shifted x-signal, i.e., the future of the x-signal by one sample.
	//%
	//% Y is a cell-array of arbitrary dimension which contains q trials of 
	//% the already embedded Y-signal.
	//%
	//% Z is a 2-dimensional p x q cell-array which contains q trials 
	//% of each of the p already embedded Z-signals. 
	//% Default empty cell array.
	//%
	//% SIGMA determines the width of the time window in samples which is
	//% used for the estimation of multivariate transfer entropy at each time 
	//% instant. Larger windows give smaller errors, but less sensitivity
	//% to temporal changes in transfer entropy, and vice versa.
	//% Default 5.
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
	//% The cell arrays X and Y are treated as 1-dimensional arrays.
	//% Each signal is a real (m x n)-matrix that contains n samples of an
	//% m-dimensional signal. The signals need not have the same dimension,
	//% but the dimension must be the same among the trials of a given signal.
	//% If the number of samples varies with each signal, the function uses 
	//% the minimum sample count among the signals. The number of trials q 
	//% must coincide with all the signals in X, Y, and Z.

	// X-signals

	const integer xSignals = mxGetNumberOfElements(inputSet[0]);
	std::vector<SignalPtr> xEnsemble;
	xEnsemble.reserve(xSignals);

	for (integer i = 0;i < xSignals;++i)
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

	// W-signals

	const integer wSignals = mxGetNumberOfElements(inputSet[1]);
	std::vector<SignalPtr> wEnsemble;
	wEnsemble.reserve(wSignals);

	for (integer i = 0;i < wSignals;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[1], i);

		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.
		const mwSize samples = mxGetN(signalArray);
		const mwSize dimension = mxGetM(signalArray);

		real* rawData = mxGetPr(signalArray);

		wEnsemble.push_back(SignalPtr(
			new Signal(samples, dimension, rawData)));
	}

	// Y-signals

	const integer ySignals = mxGetNumberOfElements(inputSet[2]);
	std::vector<SignalPtr> yEnsemble;
	yEnsemble.reserve(ySignals);

	for (integer i = 0;i < ySignals;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[2], i);

		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.
		const mwSize samples = mxGetN(signalArray);
		const mwSize dimension = mxGetM(signalArray);

		real* rawData = mxGetPr(signalArray);

		yEnsemble.push_back(SignalPtr(
			new Signal(samples, dimension, rawData)));
	}

	// Z-signals

	// Remember that Matlab uses column-major storing 
	// convention, while C++ uses row-major.

	const mwSize zSignals = mxGetDimensions(inputSet[3])[0];
	const mwSize trials = mxGetDimensions(inputSet[3])[1];

	Array<SignalPtr, 2> zEnsembleSet(trials, zSignals);

	integer index = 0;

	for (integer i = 0;i < trials;++i)
	{
		for (integer j = 0;j < zSignals;++j)
		{
			mxArray* signalArray = mxGetCell(inputSet[3], index);

			// It is intentional to assign the width
			// and height the wrong way. The reason
			// is the difference in storage conventions.
			const mwSize samples = mxGetN(signalArray);
			const mwSize dimension = mxGetM(signalArray);

			real* rawData = mxGetPr(signalArray);

			zEnsembleSet(i, j) = SignalPtr(
				new Signal(samples, dimension, rawData));

			++index;
		}
	}

	// Other parameters

	const integer sigma = *mxGetPr(inputSet[4]);
	const integer kNearest = *mxGetPr(inputSet[5]);
	const integer threads = *mxGetPr(inputSet[6]);

	omp_set_num_threads(threads);

	std::vector<real> estimateSet;

	/*
	transferEntropy(
		xEnsemble, wEnsemble, yEnsemble, 
		zEnsembleSet, sigma, kNearest, estimateSet);
	*/

	const integer samples = estimateSet.size();

	outputSet[0] = mxCreateDoubleMatrix(1, samples, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);
	
	std::copy(
		estimateSet.begin(), estimateSet.end(),
		rawResult);
}
