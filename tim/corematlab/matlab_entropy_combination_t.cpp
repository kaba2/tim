// Description: entropy_combination_t
// DocumentationOf: entropy_combination_t.m

#include "tim/corematlab/tim_matlab.h"

#include "tim/core/entropy_combination_t.h"

void force_linking_entropy_combination_t() {};

using namespace Tim;

namespace
{

	void matlabTemporalEntropyCombination(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			SignalSet,
			RangeSet,
			TimeWindowRadius,
			LagSet,
			KNearest,
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

		Array<MatlabMatrix<dreal>> signalSet = matlabAsMatrixArray<dreal>(inputSet[SignalSet]);
		MatlabMatrix<integer> lagSet = matlabAsVectorizedMatrix<integer>(inputSet[LagSet]);
		integer timeWindowRadius = matlabAsScalar<integer>(inputSet[TimeWindowRadius]);
		integer kNearest = matlabAsScalar<integer>(inputSet[KNearest]);
		MatlabMatrix<dreal> filter = matlabAsVectorizedMatrix<dreal>(inputSet[FilterIndex]);
		MatlabMatrix<dreal> rangeArray = matlabAsVectorizedMatrix<dreal>(inputSet[RangeSet]);

		integer marginals = rangeArray.rows();

		std::vector<Integer3> rangeSet;
		rangeSet.reserve(marginals);
		{
			for (integer i = 0;i < marginals;++i)
			{
				// FIX: Real weights should not be rounded to integers.

				// On Matlab's side, the range is given in the form [a, b].
				// This is the same as the range [a, b + 1[. However,
				// since Matlab indices are 1-based, this finally comes out
				// as [a - 1, b[.
				rangeSet.push_back(
					Integer3(
					rangeArray.view()(i, 0) - 1,
					rangeArray.view()(i, 1),
					rangeArray.view()(i, 2)));
			}
		}

		Signal estimate = temporalEntropyCombination(
			asSignalArray(signalSet),
			rangeSet,
			timeWindowRadius,
			lagSet.view().span(),
			kNearest,
			filter.view().span());

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
			"entropy_combination_t",
			matlabTemporalEntropyCombination);
	}

	CallFunction run(addFunction);

}
