// Description: differential_entropy_sp
// DocumentationOf: differential_entropy_sp.m

#include "tim/corematlab/tim_matlab.h"

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

		std::vector<MatlabMatrix<dreal>> xMatrices = matlabAsMatrixRange<dreal>(inputSet[X]) | ranges::to_vector;
		std::vector<Signal> xSignals = matlabMatricesAsSignals(xMatrices) | ranges::to_vector;

		dreal* outResult = matlabCreateScalar<dreal>(outputSet[Estimate]);
		*outResult = differentialEntropySp(xSignals);
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_sp",
			matlabDifferentialEntropySp);
	}

	CallFunction run(addFunction);

}
