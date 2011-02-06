// Description: differential_entropy_sp
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/differential_entropy_sp.h"

void force_linking_differential_entropy_sp() {};

using namespace Tim;

namespace
{

	void matlabDifferentialEntropySp(
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
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		std::vector<SignalPtr> xEnsemble;
		getSignals(inputSet[X], std::back_inserter(xEnsemble));

		real* outResult = createScalar<real>(outputSet[Estimate]);
		*outResult = differentialEntropySp(
			range(xEnsemble.begin(), xEnsemble.end()));
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_sp",
			matlabDifferentialEntropySp);
	}

	CallFunction run(addFunction);

}
