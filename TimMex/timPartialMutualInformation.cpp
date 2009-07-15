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
	//% A partial mutual information estimate I(X, Y | Z) from samples.
	//%
	//% I = partial_mutual_information(X, Y, Z, k)
	//%
	//% where
	//%
	//% X is a real (m x n)-matrix that contains n samples of an
	//% m-dimensional signal.
	//%
	//% Y is a real (p x n)-matrix that contains n samples of a
	//% p-dimensional signal.
	//%
	//% Z is a real (q x n)-matrix that contains n samples of a
	//% q-dimensional signal.
	//%
	//% K determines which k:th nearest neighbor the algorithm
	//% uses for estimation. Default 1.

	omp_set_num_threads(4);

	// It is intentional to assign the width
	// and height the wrong way. The reason
	// is that Matlab uses column-major storage
	// while we use row-major storage.
	const mwSize samples = mxGetN(inputSet[0]);

	const mwSize xDimension = mxGetM(inputSet[0]);
	real* xRawData = mxGetPr(inputSet[0]);

	const mwSize yDimension = mxGetM(inputSet[1]);
	real* yRawData = mxGetPr(inputSet[1]);

	const mwSize zDimension = mxGetM(inputSet[2]);
	real* zRawData = mxGetPr(inputSet[2]);
	
	integer kNearest = *mxGetPr(inputSet[3]);

	const SignalPtr xSignal = SignalPtr(
		new Signal(samples, xDimension, xRawData));
	
	const SignalPtr ySignal = SignalPtr(
		new Signal(samples, yDimension, yRawData));

	const SignalPtr zSignal = SignalPtr(
		new Signal(samples, zDimension, zRawData));

	outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	*rawResult = partialMutualInformation(
		xSignal, ySignal, zSignal, kNearest);
}
