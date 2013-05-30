// Description: entropy_combination_t
// DocumentationOf: entropy_combination_t.m

#include "tim/matlab/tim_matlab.h"

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

		Array<SignalPtr> signalSet;
		getSignalArray(inputSet[SignalSet], signalSet);

		std::vector<integer> lagSet;
		getScalars(inputSet[LagSet], std::back_inserter(lagSet));

		Array<real> rangeArray =
			asArray<real>(inputSet[RangeSet]);

		const integer marginals = rangeArray.height();

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
					rangeArray(0, i) - 1,
					rangeArray(1, i),
					rangeArray(2, i)));
			}
		}

		const integer timeWindowRadius = asScalar<integer>(inputSet[TimeWindowRadius]);
		const integer kNearest = asScalar<integer>(inputSet[KNearest]);

		std::vector<real> filter;
		getScalars(inputSet[FilterIndex], std::back_inserter(filter));

		const SignalPtr estimate = temporalEntropyCombination(
			signalSet,
			range(rangeSet.begin(), rangeSet.end()),
			timeWindowRadius,
			range(lagSet.begin(), lagSet.end()),
			kNearest,
			range(filter.begin(), filter.end()));

		const integer nans = std::max(estimate->t(), (integer)0);
		const integer skip = std::max(-estimate->t(), (integer)0); 
		const integer samples = std::max(nans + estimate->samples() - skip, (integer)0);

		Array<real> result = createArray<real>(
			Vector2i(samples, 1), outputSet[Estimate]);
		std::fill_n(result.begin(), nans, nan<real>());
		std::copy(estimate->data().begin() + skip, 
			estimate->data().end(), result.begin() + nans);
	}

	void addFunction()
	{
		matlabAddFunction(
			"entropy_combination_t",
			matlabTemporalEntropyCombination);
	}

	CallFunction run(addFunction);

}
