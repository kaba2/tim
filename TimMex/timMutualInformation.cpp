#include "mex.h"

#include "tim/core/mutual_information.h"

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

	/*
	TIMMUTUALINFORMATION A mutual information estimate from samples.
	H = timMutualInformation(D, epsilon, k)
	where
	'D' is a real (m x n)-matrix that contains n samples of an
	m-dimensional probability distribution.
	'epsilon' is the maximum relative error in distance that
	nearest neighbor searching is allowed to result in.
	Higher tolerances result in enhanced performance, but
	increases errors in the estimate. 'epsilon' defaults to 0.
	'k' determines which k:th nearest neighbor the algorithm
	uses for estimation. 'k' defaults to 1.
	*/

	// Check parameters

	if (inputs > 3) 
	{
		mexErrMsgTxt("Too many input parameters");
	} 

	if (inputs < 1) 
	{
		mexErrMsgTxt("Not enough input parameters");
	} 

	if (outputs > 1) 
	{
		mexErrMsgTxt("Too many output arguments.");
	}

	if (!mxIsDouble(inputSet[0]) || mxIsComplex(inputSet[0])) 
	{
		mexErrMsgTxt("'D' must be real and non-complex.");
	}

	if (inputs > 1)
	{
		if (mxGetM(inputSet[1]) != 1 ||
			mxGetN(inputSet[1]) != 1 ||
			!mxIsDouble(inputSet[1]) ||
			mxIsComplex(inputSet[1]))
		{
			mexErrMsgTxt("'epsilon' must be a scalar.");
		}
	}

	if (inputs > 2)
	{
		if (mxGetM(inputSet[2]) != 1 ||
			mxGetN(inputSet[2]) != 1 ||
			!mxIsDouble(inputSet[2]) ||
			mxIsComplex(inputSet[2]))
		{
			mexErrMsgTxt("'k' must be a scalar integer.");
		}
	}

	const mwSize width = mxGetM(inputSet[0]);
	const mwSize height = mxGetN(inputSet[0]);

	real* rawData = mxGetPr(inputSet[0]);
	real maxRelativeError = 0;
	if (inputs > 1)
	{
		maxRelativeError = *mxGetPr(inputSet[1]);
		if (maxRelativeError < 0)
		{
			mexErrMsgTxt("'epsilon' must be positive'.");
		}
	}

	integer kNearest = 1;
	if (inputs > 2)
	{
		kNearest = *mxGetPr(inputSet[2]);
		if (kNearest < 1)
		{
			mexErrMsgTxt("'k' must be at least 1.");
		}
	}

	const SignalPtr data = SignalPtr(
		new Signal(height, width, withAliasing(rawData)));

	/*
	printf("k = %d, epsilon = %f\n", kNearest, maxRelativeError);
	for (integer y = 0;y < height;++y)
	{
		for (integer x = 0;x < width;++x)
		{
			printf("%f ", *(rawData + y * width + x));
		}
		printf("\n");
	}
	*/
	
	outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	*rawResult = mutualInformation(data, kNearest, 
		maxRelativeError, InfinityNormBijection<real>());
}
