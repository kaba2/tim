// Description: tsallis_entropy_lps
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/tsallis_entropy_lps.h"

using namespace Tim;

namespace
{

	void matlabTsallisEntropyLps(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum
		{
			xIndex,
			qIndex,
			kNearestSuggestionIndex,
			threadsIndex
		};

		std::vector<SignalPtr> xEnsemble;
		getSignals(inputSet[xIndex], std::back_inserter(xEnsemble));

		const real q = asReal(inputSet[qIndex]);
		const integer kNearestSuggestion = asInteger(inputSet[kNearestSuggestionIndex]);
		const integer threads = asInteger(inputSet[threadsIndex]);
		setNumberOfThreads(threads);

		outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
		real* rawResult = mxGetPr(outputSet[0]);

		*rawResult = tsallisEntropyLps(
			range(xEnsemble.begin(), xEnsemble.end()), 
			q, kNearestSuggestion);
	}

	void addFunction()
	{
		matlabAddFunction(
			"tsallis_entropy_lps",
			matlabTsallisEntropyLps);
	}

	CallFunction run(addFunction);

}
