// Description: mutual_information_naive
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/mutual_information_naive.h"

void force_linking_mutual_information_naive() {};

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

		const SignalPtr data = asSignal(inputSet[xIndex]);
		const integer bins = asInteger(inputSet[binsIndex]);

		const integer dimension = data->dimension();

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

