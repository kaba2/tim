// Description: differential_entropy_nk
// DocumentationOf: differential_entropy_nk.m

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
		ENSURE_OP(outputs, <=, Outputs);

		std::vector<Signal> xEnsemble = getSignals(inputSet[X]);

		integer intrinsicDimension = 0;

		real entropy = 
			differentialEntropyNk(
			countingRange(xEnsemble.begin(), xEnsemble.end()), 
			Euclidean_NormBijection<real>(),

			&intrinsicDimension);

		if (outputs > 0)
		{
			real* outEntropy =
				createScalar<real>(outputSet[Estimate]);
			*outEntropy = entropy;
		}

		if (outputs > 1)
		{
			integer* outIntrinsicDimension =
				createScalar<integer>(outputSet[IntrinsicDimension]);
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
