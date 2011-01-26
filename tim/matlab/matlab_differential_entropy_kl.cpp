// Description: differential_entropy_kl
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/differential_entropy_kl.h"

void force_linking_differential_entropy_kl() {};

using namespace Tim;

namespace
{

	void matlabDifferentialEntropyKl(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			X,
			KNearest,
			Inputs
		};

		enum Output
		{
			Estimate,
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		std::vector<SignalPtr> xEnsemble;
		getSignals(inputSet[X], std::back_inserter(xEnsemble));

		const integer kNearest = asInteger(inputSet[KNearest]);

		outputSet[Estimate] = mxCreateDoubleMatrix(1, 1, mxREAL);
		real* rawResult = mxGetPr(outputSet[Estimate]);

		*rawResult = differentialEntropyKl(
			range(xEnsemble.begin(), xEnsemble.end()), 
			kNearest);
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_kl",
			matlabDifferentialEntropyKl);
	}

	CallFunction run(addFunction);

}
