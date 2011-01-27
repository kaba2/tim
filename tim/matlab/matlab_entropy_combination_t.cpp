// Description: entropy_combination_t
// Documentation: tim_matlab_functions.txt

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

		const integer marginals = mxGetM(inputSet[RangeSet]);
		std::vector<Integer3> rangeSet;
		rangeSet.reserve(marginals);
		{
			real* rawData = mxGetPr(inputSet[RangeSet]);
			for (integer i = 0;i < marginals;++i)
			{
				rangeSet.push_back(Integer3(*rawData - 1, *(rawData + marginals), *(rawData + 2 * marginals)));
				//printf("%d %d %d ", rangeSet.back()[0], rangeSet.back()[1], rangeSet.back()[1]);
				++rawData;
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

		const integer nans = std::max(estimate->t(), 0);
		const integer skip = std::max(-estimate->t(), 0); 
		const integer samples = std::max(nans + estimate->samples() - skip, 0);

		outputSet[Estimate] = mxCreateDoubleMatrix(1, samples, mxREAL);
		real* rawResult = mxGetPr(outputSet[Estimate]);
		std::fill_n(rawResult, nans, nan<real>());
		std::copy(estimate->data().begin() + skip, 
			estimate->data().end(), rawResult + nans);
	}

	void addFunction()
	{
		matlabAddFunction(
			"entropy_combination_t",
			matlabTemporalEntropyCombination);
	}

	CallFunction run(addFunction);

}
