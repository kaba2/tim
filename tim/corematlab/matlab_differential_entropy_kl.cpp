// Description: differential_entropy_kl
// DocumentationOf: differential_entropy_kl.m

#include "tim/corematlab/tim_matlab.h"

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

		std::vector<MatlabMatrix<dreal>> xMatrices = matlabAsMatrixRange<dreal>(inputSet[X]) | ranges::to_vector;
		std::vector<Signal> xSignals = matlabMatricesAsSignals(xMatrices) | ranges::to_vector;
		
		integer kNearest = matlabAsScalar<integer>(inputSet[KNearest]);

		dreal* outResult = matlabCreateScalar<dreal>(outputSet[Estimate]);
		*outResult = differentialEntropyKl(
			xSignals,
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
