// Description: Differential entropy of a generalized normal distribution
// DocumentationOf: differential_entropy_uniform.m

#include "tim/matlab/tim_matlab.h"

#include "tim/core/differential_entropy_generalized_normal.h"

void force_linking_differential_entropy_generalized_normal() {};

using namespace Tim;

namespace
{

	void matlabDifferentialEntropyGeneralizedNormal(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			Dimension,
			Shape,
			Scale,
			Inputs
		};

		enum Output
		{
			Entropy,
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		integer dimension =
			asScalar<integer>(inputSet[Dimension]);
		real shape = 
			asScalar<real>(inputSet[Shape]);
		real scale = 
			asScalar<real>(inputSet[Scale]);

		real* entropy = createScalar<real>(outputSet[Entropy]);
		*entropy = differentialEntropyGeneralizedNormal<real>(dimension, shape, scale);
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_generalized_normal",
			matlabDifferentialEntropyGeneralizedNormal);
	}

	CallFunction run(addFunction);

}
