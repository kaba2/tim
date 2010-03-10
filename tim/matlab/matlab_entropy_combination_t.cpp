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

		Array<SignalPtr, 2> signalSet;
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

		const integer outputWidth = minSamples(
			forwardRange(signalSet.begin(), signalSet.end()));

		outputSet[0] = mxCreateDoubleMatrix(1, outputWidth, mxREAL);
		real* rawResult = mxGetPr(outputSet[0]);

		temporalEntropyCombination(
			signalSet,
			forwardRange(rangeSet.begin(), rangeSet.end()),
			timeWindowRadius,
			rawResult,
			forwardRange(lagSet.begin(), lagSet.end()),
			kNearest,
			forwardRange(filter.begin(), filter.end()));
	}

	void addFunction()
	{
		matlabAddFunction(
			"entropy_combination_t",
			matlabTemporalEntropyCombination);
	}

	CallFunction run(addFunction);

}
