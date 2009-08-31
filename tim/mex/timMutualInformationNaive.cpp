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

	enum
	{
		xIndex,
		binsIndex
	};

	const integer samples = mxGetN(inputSet[xIndex]);
	const integer dimension = mxGetM(inputSet[xIndex]);

	real* rawData = mxGetPr(inputSet[xIndex]);
	integer bins = *mxGetPr(inputSet[binsIndex]);

	const SignalPtr data = SignalPtr(
		new Signal(samples, dimension, withAliasing(rawData)));
	
	outputSet[0] = mxCreateDoubleMatrix(dimension, dimension, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	MatrixD result(dimension, dimension, withAliasing(rawResult));

	mutualInformationFromBinning(data, bins, result);
}
