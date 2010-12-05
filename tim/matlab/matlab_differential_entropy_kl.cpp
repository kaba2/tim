// Description: differential_entropy_kl
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/differential_entropy_kl.h"

using namespace Tim;

namespace
{

	void matlabDifferentialEntropyKl(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum
		{
			xIndex,
			kNearestIndex,
			threadsIndex
		};

		std::vector<SignalPtr> xEnsemble;
		getSignals(inputSet[xIndex], std::back_inserter(xEnsemble));

		const integer kNearest = getInteger(inputSet[kNearestIndex]);
		const integer threads = getInteger(inputSet[threadsIndex]);
		setNumberOfThreads(threads);

		outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
		real* rawResult = mxGetPr(outputSet[0]);

		*rawResult = differentialEntropyKl(
			range(xEnsemble.begin(), xEnsemble.end()), 
			kNearest);
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_kl",
			matlabDifferentialEntropyKl);
	}

	CallFunction run(addFunction);

}
