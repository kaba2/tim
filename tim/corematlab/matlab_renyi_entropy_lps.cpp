// Description: renyi_entropy_lps
// DocumentationOf: renyi_entropy_lps.m

#include "tim/corematlab/tim_matlab.h"

#include "tim/core/renyi_entropy_lps.h"

void force_linking_renyi_entropy_lps() {};

using namespace Tim;

namespace
{

	void matlabRenyiEntropyLps(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			X,
			Q,
			KNearestSuggestion,
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

		dreal q = matlabAsScalar<dreal>(inputSet[Q]);
		integer kNearestSuggestion = matlabAsScalar<integer>(inputSet[KNearestSuggestion]);


		dreal* outResult = matlabCreateScalar<dreal>(outputSet[Estimate]);
		*outResult = renyiEntropyLps(
			xSignals,
			q, kNearestSuggestion);
	}

	void addFunction()
	{
		matlabAddFunction(
			"renyi_entropy_lps",
			matlabRenyiEntropyLps);
	}

	CallFunction run(addFunction);

}
