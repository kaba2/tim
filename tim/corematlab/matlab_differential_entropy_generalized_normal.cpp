// Description: Differential entropy of a generalized normal distribution
// DocumentationOf: differential_entropy_uniform.m

#include "tim/corematlab/tim_matlab.h"

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
			matlabAsScalar<integer>(inputSet[Dimension]);
		dreal shape = 
			matlabAsScalar<dreal>(inputSet[Shape]);
		dreal scale = 
			matlabAsScalar<dreal>(inputSet[Scale]);

		dreal* entropy = matlabCreateScalar<dreal>(outputSet[Entropy]);
		*entropy = differentialEntropyGeneralizedNormal<dreal>(dimension, shape, scale);
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_generalized_normal",
			matlabDifferentialEntropyGeneralizedNormal);
	}

	CallFunction run(addFunction);

}
