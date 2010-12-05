// Description: entropy_combination
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/entropy_combination.h"

using namespace Tim;

namespace
{

	void matlabEntropyCombination(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum
		{
			signalSetIndex,
			rangeSetIndex,
			lagSetIndex,
			kNearestIndex,
			threadsIndex
		};

		Array<SignalPtr> signalSet;
		getSignalArray(inputSet[signalSetIndex], signalSet);

		std::vector<integer> lagSet;
		getIntegers(inputSet[lagSetIndex], std::back_inserter(lagSet));

		const integer marginals = mxGetM(inputSet[rangeSetIndex]);
		//printf("%d marginals", marginals);
		std::vector<Integer3> rangeSet;
		rangeSet.reserve(marginals);
		{
			real* rawData = mxGetPr(inputSet[rangeSetIndex]);
			for (integer i = 0;i < marginals;++i)
			{
				rangeSet.push_back(Integer3(*rawData - 1, *(rawData + marginals), *(rawData + 2 * marginals)));
				//printf("%d %d %d ", rangeSet.back()[0], rangeSet.back()[1], rangeSet.back()[2]);
				++rawData;
			}
		}

		const integer kNearest = getInteger(inputSet[kNearestIndex]);
		const integer threads = getInteger(inputSet[threadsIndex]);
		setNumberOfThreads(threads);

		outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
		real* rawResult = mxGetPr(outputSet[0]);

		*rawResult = entropyCombination(
			signalSet,
			range(rangeSet.begin(), rangeSet.end()),
			range(lagSet.begin(), lagSet.end()),
			kNearest);
	}

	void addFunction()
	{
		matlabAddFunction(
			"entropy_combination",
			matlabEntropyCombination);
	}

	CallFunction run(addFunction);

}
