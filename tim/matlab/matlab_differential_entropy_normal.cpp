// Description: Differential entropy of a normal distribution
// DocumentationOf: differential_entropy_normal.m

#include "tim/matlab/tim_matlab.h"

#include "tim/core/differential_entropy_normal.h"

void force_linking_differential_entropy_normal() {};

using namespace Tim;

namespace
{

	void matlabDifferentialEntropyNormal(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			Dimension,
			CovarianceDeterminant,
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
		real covarianceDeterminant = 
			asScalar<real>(inputSet[CovarianceDeterminant]);

		real* entropy = 
			createScalar<real>(outputSet[Entropy]);
		*entropy = 
			differentialEntropyNormal<real>(dimension, covarianceDeterminant);
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_normal",
			matlabDifferentialEntropyNormal);
	}

	CallFunction run(addFunction);

}
