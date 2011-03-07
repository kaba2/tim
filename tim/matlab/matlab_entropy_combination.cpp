// Description: entropy_combination
// DocumentationOf: entropy_combination.m

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
		getScalars(inputSet[LagSet], std::back_inserter(lagSet));

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

		const integer kNearest = asScalar<integer>(inputSet[KNearest]);

		real* outResult = createScalar<real>(outputSet[Estimate]);
		*outResult = entropyCombination(
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
