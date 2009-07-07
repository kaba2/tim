#include "mex.h"

#include "tim/core/partial_mutual_information.h"

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

	// It is intentional to assign the width
	// and height the wrong way. The reason
	// is that Matlab uses column-major storage
	// while we use row-major storage.
	const mwSize samples = mxGetN(inputSet[0]);

	const mwSize xDimension = mxGetN(inputSet[0]);
	real* xRawData = mxGetPr(inputSet[0]);

	const mwSize yDimension = mxGetN(inputSet[1]);
	real* yRawData = mxGetPr(inputSet[1]);

	const mwSize zDimension = mxGetN(inputSet[2]);
	real* zRawData = mxGetPr(inputSet[2]);
	
	integer kNearest = *mxGetPr(inputSet[3]);

	const SignalPtr xSignal = SignalPtr(
		new Signal(xDimension, samples, xRawData));
	
	const SignalPtr ySignal = SignalPtr(
		new Signal(yDimension, samples, yRawData));

	const SignalPtr zSignal = SignalPtr(
		new Signal(zDimension, samples, zRawData));

	outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	*rawResult = partialMutualInformation(
		xSignal, ySignal, zSignal, kNearest);
}
