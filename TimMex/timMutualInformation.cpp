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

	// It is intentional to assign the width
	// and height the wrong way. The reason
	// is that Matlab uses column-major storage
	// while we use row-major storage.
	const mwSize width = mxGetM(inputSet[0]);
	const mwSize height = mxGetN(inputSet[0]);

	real* rawData = mxGetPr(inputSet[0]);
	real maxRelativeError = *mxGetPr(inputSet[1]);
	integer kNearest = *mxGetPr(inputSet[2]);

	const SignalPtr data = SignalPtr(
		new Signal(height, width, rawData));
	
	outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	*rawResult = mutualInformation(data, kNearest, 
		maxRelativeError, InfinityNormBijection<real>());
}
