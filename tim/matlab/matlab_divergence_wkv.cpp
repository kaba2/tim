// Description: divergence_wkv
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/divergence_wkv.h"

using namespace Tim;

namespace
{

	void matlabDivergenceWkv(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum
		{
			xIndex,
			yIndex,
			threadsIndex
		};

		std::vector<SignalPtr> xEnsemble;
		getSignals(inputSet[xIndex], std::back_inserter(xEnsemble));

		std::vector<SignalPtr> yEnsemble;
		getSignals(inputSet[yIndex], std::back_inserter(yEnsemble));

		const integer threads = asInteger(inputSet[threadsIndex]);
		setNumberOfThreads(threads);

		outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
		real* rawResult = mxGetPr(outputSet[0]);

		*rawResult = divergenceWkv(
			range(xEnsemble.begin(), xEnsemble.end()), 
			range(yEnsemble.begin(), yEnsemble.end()));
	}

	void addFunction()
	{
		matlabAddFunction(
			"divergence_wkv",
			matlabDivergenceWkv);
	}

	CallFunction run(addFunction);

}
