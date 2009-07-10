#include "mex.h"

#include "tim/core/mutual_information.h"

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

	//% MUTUAL_INFORMATION A mutual information estimate from samples.
	//% I = mutual_information(S, k, threads)
	//% where
	//% S is a cell array of arbitrary dimension that contains p signals 
	//% (it is addressed as a 1d cell-array).
	//% Each signal is a real (m x n)-matrix that contains n samples of an
	//% m-dimensional signal. If the number of samples varies with each
	//% signal, the function uses the minimum sample count among the signals.
	//% K determines which k:th nearest neighbor the algorithm
	//% uses for estimation. Default 1.
	//% THREADS determines the number of threads to use for parallelization.
	//% A good value is to set the number of threads to the number of
	//% cores a machine has. To keep the machine responsive, you might
	//% choose to spare one core for other work. Default 1 (no parallelization).

	const integer signals = mxGetNumberOfElements(inputSet[0]);
	std::vector<SignalPtr> signalSet;
	signalSet.reserve(signals);

	for (integer i = 0;i < signals;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[0], i);

		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.
		const mwSize samples = mxGetN(signalArray);
		const mwSize dimension = mxGetM(signalArray);

		real* rawData = mxGetPr(signalArray);

		signalSet.push_back(SignalPtr(
			new Signal(samples, dimension, rawData)));
	}

	const integer kNearest = *mxGetPr(inputSet[1]);
	const integer threads = *mxGetPr(inputSet[2]);

	omp_set_num_threads(threads);

	outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	*rawResult = mutualInformation(signalSet, kNearest);
}
