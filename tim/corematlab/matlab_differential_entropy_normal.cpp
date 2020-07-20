// Description: Differential entropy of a normal distribution
// DocumentationOf: differential_entropy_normal.m

#include "tim/corematlab/tim_matlab.h"

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
			matlabAsScalar<integer>(inputSet[Dimension]);
		dreal covarianceDeterminant = 
			matlabAsScalar<dreal>(inputSet[CovarianceDeterminant]);

		dreal* entropy = 
			matlabCreateScalar<dreal>(outputSet[Entropy]);
		*entropy = 
			differentialEntropyNormal<dreal>(dimension, covarianceDeterminant);
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_normal",
			matlabDifferentialEntropyNormal);
	}

	CallFunction run(addFunction);

}
