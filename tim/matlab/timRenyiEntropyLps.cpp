// Description: renyi_entropy_lps
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/renyi_entropy_lps.h"

using namespace Tim;

namespace
{

	void matlabRenyiEntropyLps(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum
		{
			xIndex,
			qIndex,
			maxRelativeErrorIndex,
			kNearestSuggestionIndex,
			threadsIndex
		};

		std::vector<SignalPtr> xEnsemble;
		getSignals(inputSet[xIndex], std::back_inserter(xEnsemble));

		const real q = getReal(inputSet[qIndex]);
		const real maxRelativeError = getReal(inputSet[maxRelativeErrorIndex]);
		const integer kNearestSuggestion = getInteger(inputSet[kNearestSuggestionIndex]);
		const integer threads = getInteger(inputSet[threadsIndex]);
		setNumberOfThreads(threads);

		outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
		real* rawResult = mxGetPr(outputSet[0]);

		*rawResult = renyiEntropyLps(
			randomAccessRange(xEnsemble.begin(), xEnsemble.end()), 
			q, maxRelativeError, kNearestSuggestion);
	}

	void addFunction()
	{
		matlabAddFunction(
			"renyi_entropy_lps",
			matlabRenyiEntropyLps);
	}

	CallFunction run(addFunction);

}
