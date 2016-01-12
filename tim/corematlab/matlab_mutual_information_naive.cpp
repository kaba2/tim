// Description: mutual_information_naive
// DocumentationOf: mutual_information_naive.m

#include "tim/corematlab/tim_matlab.h"

#include "tim/core/mutual_information_naive.h"

void force_linking_mutual_information_naive() {};

using namespace Tim;

namespace
{

	void matlabMutualInformationNaive(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			X,
			Bins,
			Inputs
		};
		
		enum Output
		{
			Estimate,
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		const Signal& data = asSignal(inputSet[X]);
		integer bins = matlabAsScalar<integer>(inputSet[Bins]);

		integer n = data.dimension();

		Array<real> result =
			matlabCreateArray<real>(Vector2i(n, n), 
			outputSet[Estimate]);

		result = mutualInformationFromBinning(data, bins);
	}

	void addFunction()
	{
		matlabAddFunction(
			"mutual_information_naive",
			matlabMutualInformationNaive);
	}

	CallFunction run(addFunction);

}

