// Description: differential_entropy_nk
// DocumentationOf: differential_entropy_nk.m

#include "tim/corematlab/tim_matlab.h"

#include "tim/core/differential_entropy_nk.h"

#include <pastel/math/normbijection/euclidean_normbijection.h>

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
		ENSURE_OP(outputs, <=, Outputs);

		std::vector<MatlabMatrix<dreal>> xMatrices = matlabAsMatrixRange<dreal>(inputSet[X]) | ranges::to_vector;
		std::vector<Signal> xSignals = matlabMatricesAsSignals(xMatrices) | ranges::to_vector;

		integer intrinsicDimension = 0;

		dreal entropy = 
			differentialEntropyNk(
			xSignals, 
			Euclidean_Norm<dreal>(),

			&intrinsicDimension);

		if (outputs > 0)
		{
			dreal* outEntropy =
				matlabCreateScalar<dreal>(outputSet[Estimate]);
			*outEntropy = entropy;
		}

		if (outputs > 1)
		{
			integer* outIntrinsicDimension =
				matlabCreateScalar<integer>(outputSet[IntrinsicDimension]);
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
