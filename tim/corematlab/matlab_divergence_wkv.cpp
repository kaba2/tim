// Description: divergence_wkv
// DocumentationOf: divergence_wkv.m

#include "tim/corematlab/tim_matlab.h"

#include "tim/core/divergence_wkv.h"

void force_linking_divergence_wkv() {};

using namespace Tim;

namespace
{

	void matlabDivergenceWkv(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			X,
			Y,
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

		std::vector<MatlabMatrix<dreal>> yMatrices = matlabAsMatrixRange<dreal>(inputSet[Y]) | ranges::to_vector;
		std::vector<Signal> ySignals = matlabMatricesAsSignals(yMatrices) | ranges::to_vector;

		dreal* outResult = matlabCreateScalar<dreal>(outputSet[Estimate]);
		*outResult = divergenceWkv(
			xSignals, 
			ySignals);
	}

	void addFunction()
	{
		matlabAddFunction(
			"divergence_wkv",
			matlabDivergenceWkv);
	}

	CallFunction run(addFunction);

}
