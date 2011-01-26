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
		enum Input
		{
			X,
			Bins,
			Inputs
		};
		
		enum Output
		{
			Estimate,
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		const SignalPtr data = asSignal(inputSet[X]);
		const integer bins = asInteger(inputSet[Bins]);

		const integer dimension = data->dimension();

		outputSet[Estimate] = 
			mxCreateDoubleMatrix(dimension, dimension, mxREAL);
		real* rawResult = mxGetPr(outputSet[Estimate]);

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

