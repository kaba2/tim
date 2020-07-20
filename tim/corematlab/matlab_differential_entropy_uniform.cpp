// Description: Differential entropy of a uniform distribution
// DocumentationOf: differential_entropy_uniform.m

#include "tim/corematlab/tim_matlab.h"

#include "tim/core/differential_entropy_uniform.h"

void force_linking_differential_entropy_uniform() {};

using namespace Tim;

namespace
{

	void matlabDifferentialEntropyUniform(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			SupportVolume,
			Inputs
		};

		enum Output
		{
			Entropy,
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		dreal supportVolume = 
			matlabAsScalar<dreal>(inputSet[SupportVolume]);

		dreal* entropy = matlabCreateScalar<dreal>(outputSet[Entropy]);
		*entropy = differentialEntropyUniform<dreal>(supportVolume);
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_uniform",
			matlabDifferentialEntropyUniform);
	}

	CallFunction run(addFunction);

}
