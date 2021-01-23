// Description: tsallis_entropy_lps_t
// DocumentationOf: tsallis_entropy_lps_t.m

#include "tim/corematlab/tim_matlab.h"

#include "tim/core/tsallis_entropy_lps.h"

void force_linking_tsallis_entropy_lps_t() {};

using namespace Tim;

namespace
{

	void matlabTemporalTsallisEntropyLps(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			X,
			TimeWindowRadius,
			Q,
			KNearestSuggestion,
			FilterIndex,
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

		integer timeWindowRadius = matlabAsScalar<integer>(inputSet[TimeWindowRadius]);
		dreal q = matlabAsScalar<dreal>(inputSet[Q]);
		integer kNearestSuggestion = matlabAsScalar<integer>(inputSet[KNearestSuggestion]);

		std::vector<dreal> filter;
		matlabGetScalars(inputSet[FilterIndex], std::back_inserter(filter));

		SignalData estimate = temporalTsallisEntropyLps(
			xSignals,
			timeWindowRadius, 
			q,
			kNearestSuggestion,
			range(std::begin(filter), std::end(filter)));

		integer nans = std::max(estimate.t(), (integer)0);
		integer skip = std::max(-estimate.t(), (integer)0); 
		integer samples = std::max(nans + estimate.samples() - skip, (integer)0);

		MatrixView<dreal> result = matlabCreateMatrix<dreal>(1, samples, outputSet[Estimate]);
		ranges::fill(result.slicex(0, nans).range(), (dreal)Nan());
		ranges::copy(
			estimate.data().slicex(skip).range(),
			std::begin(result.slicex(nans).range()));

	}

	void addFunction()
	{
		matlabAddFunction(
			"tsallis_entropy_lps_t",
			matlabTemporalTsallisEntropyLps);
	}

	CallFunction run(addFunction);

}
