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
		printf("%d x %d signal-set.", signalSet.width(), signalSet.height());

		std::vector<integer> lagSet;
		getScalars(inputSet[LagSet], std::back_inserter(lagSet));
		log() << "Lags: ";
		for (integer i = 0;i < lagSet.size();++i)
		{
			log() << lagSet[i] << ", ";
		}
		log() << logNewLine;

		RealArrayPtr rangeArray =
			asArray<real>(inputSet[RangeSet]);

		const integer marginals = rangeArray->height();
		log() << marginals << " marginals." << logNewLine;

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
					(*rangeArray)(0, i) - 1,
					(*rangeArray)(1, i),
					(*rangeArray)(2, i)));
					
				printf("%d %d %d ", rangeSet.back()[0], rangeSet.back()[1], rangeSet.back()[2]);
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
