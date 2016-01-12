// Description: entropy_combination
// DocumentationOf: entropy_combination.m

#include "tim/corematlab/tim_matlab.h"

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
			Estimator,
			Inputs
		};

		enum Output
		{
			Estimate,
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		Array<Signal> signalSet = getSignalArray(inputSet[SignalSet]);

		std::vector<integer> lagSet;
		matlabGetScalars(inputSet[LagSet], std::back_inserter(lagSet));

		Array<real> rangeArray =
			matlabAsArray<real>(inputSet[RangeSet]);

		integer marginals = rangeArray.height();

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

		integer kNearest = matlabAsScalar<integer>(inputSet[KNearest]);
		std::string estimator = matlabAsString(inputSet[Estimator]);


		real* outResult = matlabCreateScalar<real>(outputSet[Estimate]);

		if (estimator == "log")
		{
			*outResult = entropyCombination(
				signalSet,
				range(rangeSet.cbegin(), rangeSet.cend()),
				range(lagSet.cbegin(), lagSet.cend()),
				kNearest,
				Log_LocalEstimator());
		}
		else if (estimator == "digamma")
		{
			*outResult = entropyCombination(
				signalSet,
				range(rangeSet.cbegin(), rangeSet.cend()),
				range(lagSet.cbegin(), lagSet.cend()),
				kNearest,
				Digamma_LocalEstimator());
		}
		else if (estimator == "digamma_density")
		{
			*outResult = entropyCombination(
				signalSet,
				range(rangeSet.cbegin(), rangeSet.cend()),
				range(lagSet.cbegin(), lagSet.cend()),
				kNearest,
				DigammaDensity_LocalEstimator());
		}
	}

	void addFunction()
	{
		matlabAddFunction(
			"entropy_combination",
			matlabEntropyCombination);
	}

	CallFunction run(addFunction);

}
