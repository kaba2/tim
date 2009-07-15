#include "mex.h"

#include "tim/core/mutual_information_naive.h"

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

	//% MUTUAL_INFORMATION_NAIVE
	//% A mutual information estimate from samples.
	//%
	//% I = mutual_information_naive(S, bins)
	//%
	//% where
	//%
	//% S is a real (m x n)-matrix that contains n samples of an
	//% m-dimensional signal.
	//%
	//% BINS determines the number of bins to use for 1d
	//% distribution estimation. Default 100.

	// It is intentional to assign the width
	// and height the wrong way. The reason
	// is that Matlab uses column-major storage
	// while we use row-major storage.
	const mwSize samples = mxGetN(inputSet[0]);
	const mwSize dimension = mxGetM(inputSet[0]);

	real* rawData = mxGetPr(inputSet[0]);
	integer bins = *mxGetPr(inputSet[1]);

	const SignalPtr data = SignalPtr(
		new Signal(samples, dimension, withAliasing(rawData)));
	
	outputSet[0] = mxCreateDoubleMatrix(dimension, dimension, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	MatrixD result(dimension, dimension, withAliasing(rawResult));

	mutualInformationFromBinning(data, bins, result);
}
