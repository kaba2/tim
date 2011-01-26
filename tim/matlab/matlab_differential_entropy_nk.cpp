// Description: differential_entropy_nk
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/differential_entropy_nk.h"

#include <pastel/math/euclidean_normbijection.h>

void force_linking_differential_entropy_nk() {};

using namespace Tim;

namespace
{

	void matlabDifferentialEntropyNk(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			X,
			Inputs
		};

		enum Output
		{
			Estimate,
			IntrinsicDimension,
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		std::vector<SignalPtr> xEnsemble;
		getSignals(inputSet[X], std::back_inserter(xEnsemble));

		integer intrinsicDimension = 0;

		const real entropy = 
			differentialEntropyNk(
			range(xEnsemble.begin(), xEnsemble.end()), 
			Euclidean_NormBijection<real>(),
			&intrinsicDimension);

		if (outputs > 0)
		{
			outputSet[Estimate] = mxCreateDoubleMatrix(1, 1, mxREAL);
			real* outEntropy = mxGetPr(outputSet[Estimate]);
			*outEntropy = entropy;
		}

		if (outputs > 1)
		{
			outputSet[IntrinsicDimension] = 
				mxCreateDoubleMatrix(1, 1, mxREAL);
			real* outIntrinsicDimension = mxGetPr(outputSet[IntrinsicDimension]);
			*outIntrinsicDimension = intrinsicDimension;
		}
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_nk",
			matlabDifferentialEntropyNk);
	}

	CallFunction run(addFunction);

}
