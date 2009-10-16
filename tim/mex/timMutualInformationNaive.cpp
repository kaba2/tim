#include "tim/mex/tim_mex.h"

#include "tim/core/mutual_information_naive.h"

using namespace Tim;

void mexFunction(int outputs, mxArray *outputSet[],
				 int inputs, const mxArray *inputSet[])
{
	enum
	{
		xIndex,
		binsIndex
	};

	const integer samples = mxGetN(inputSet[xIndex]);
	const integer dimension = mxGetM(inputSet[xIndex]);
	real* rawData = mxGetPr(inputSet[xIndex]);

	integer bins = getInteger(inputSet[binsIndex]);

	const SignalPtr data = SignalPtr(
		new Signal(samples, dimension, withAliasing(rawData)));
	
	outputSet[0] = mxCreateDoubleMatrix(dimension, dimension, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	MatrixD result(dimension, dimension, withAliasing(rawResult));

	mutualInformationFromBinning(data, bins, result);
}
