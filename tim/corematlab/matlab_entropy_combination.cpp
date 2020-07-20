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
		MatlabMatrix<integer> lagSet = matlabAsMatrix<integer>(inputSet[LagSet]);
		MatlabMatrix<dreal> rangeArray = matlabAsMatrix<dreal>(inputSet[RangeSet]);
		integer kNearest = matlabAsScalar<integer>(inputSet[KNearest]);

		integer marginals = rangeArray.rows();
		ENSURE_OP(rangeArray.cols(), ==, 3);

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

		dreal result = entropyCombination(
			asSignalArray(signalSet),
			rangeSet,
			lagSet.view().range(),
			kNearest,
			Digamma_LocalEstimator());

		*matlabCreateScalar<dreal>(outputSet[Estimate]) = result;

	}

	void addFunction()
	{
		matlabAddFunction(
			"entropy_combination",
			matlabEntropyCombination);
	}

	CallFunction run(addFunction);

}
