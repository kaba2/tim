// Description: mutual_information_naive
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/mutual_information_naive.h"

using namespace Tim;

namespace
{

	void matlabMutualInformationNaive(
		int outputs, mxArray *outputSet[],
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

	void addFunction()
	{
		matlabAddFunction(
			"mutual_information_naive",
			matlabMutualInformationNaive);
	}

	CallFunction run(addFunction);

}

