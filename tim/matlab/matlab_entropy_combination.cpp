// Description: entropy_combination
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/entropy_combination.h"

void force_linking_entropy_combination() {};

using namespace Tim;

namespace
{

	void matlabEntropyCombination(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum
		{
			SignalSet,
			RangeSet,
			LagSet,
			KNearest,
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
		getIntegers(inputSet[LagSet], std::back_inserter(lagSet));

		const integer marginals = mxGetM(inputSet[RangeSet]);
		//printf("%d marginals", marginals);
		std::vector<Integer3> rangeSet;
		rangeSet.reserve(marginals);
		{
			real* rawData = mxGetPr(inputSet[RangeSet]);
			for (integer i = 0;i < marginals;++i)
			{
				rangeSet.push_back(Integer3(*rawData - 1, *(rawData + marginals), *(rawData + 2 * marginals)));
				//printf("%d %d %d ", rangeSet.back()[0], rangeSet.back()[1], rangeSet.back()[2]);
				++rawData;
			}
		}

		const integer kNearest = asInteger(inputSet[KNearest]);

		outputSet[Estimate] = mxCreateDoubleMatrix(1, 1, mxREAL);
		real* rawResult = mxGetPr(outputSet[Estimate]);

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
