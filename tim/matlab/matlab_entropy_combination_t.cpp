// Description: entropy_combination_t
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/entropy_combination_t.h"

using namespace Tim;

namespace
{

	void matlabTemporalEntropyCombination(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum
		{
			signalSetIndex,
			rangeSetIndex,
			timeWindowRadiusIndex,
			lagSetIndex,
			kNearestIndex,
			filterIndex,
			threadsIndex
		};

		Array<SignalPtr> signalSet;
		getSignalArray(inputSet[signalSetIndex], signalSet);

		std::vector<integer> lagSet;
		getIntegers(inputSet[lagSetIndex], std::back_inserter(lagSet));

		const integer marginals = mxGetM(inputSet[rangeSetIndex]);
		std::vector<Integer3> rangeSet;
		rangeSet.reserve(marginals);
		{
			real* rawData = mxGetPr(inputSet[rangeSetIndex]);
			for (integer i = 0;i < marginals;++i)
			{
				rangeSet.push_back(Integer3(*rawData - 1, *(rawData + marginals), *(rawData + 2 * marginals)));
				//printf("%d %d %d ", rangeSet.back()[0], rangeSet.back()[1], rangeSet.back()[1]);
				++rawData;
			}
		}

		const integer timeWindowRadius = getInteger(inputSet[timeWindowRadiusIndex]);
		const integer kNearest = getInteger(inputSet[kNearestIndex]);

		std::vector<real> filter;
		getReals(inputSet[filterIndex], std::back_inserter(filter));

		const integer threads = getInteger(inputSet[threadsIndex]);
		setNumberOfThreads(threads);

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

		outputSet[0] = mxCreateDoubleMatrix(1, samples, mxREAL);
		real* rawResult = mxGetPr(outputSet[0]);
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
