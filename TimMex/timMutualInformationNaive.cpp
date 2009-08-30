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

	const integer samples = mxGetN(inputSet[0]);
	const integer dimension = mxGetM(inputSet[0]);

	real* rawData = mxGetPr(inputSet[0]);
	integer bins = *mxGetPr(inputSet[1]);

	const SignalPtr data = SignalPtr(
		new Signal(samples, dimension, withAliasing(rawData)));
	
	outputSet[0] = mxCreateDoubleMatrix(dimension, dimension, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	MatrixD result(dimension, dimension, withAliasing(rawResult));

	mutualInformationFromBinning(data, bins, result);
}
