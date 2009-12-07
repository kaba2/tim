#include "tim/mex/tim_mex.h"

#include "tim/core/differential_entropy_kl.h"

using namespace Tim;

void mexFunction(int outputs, mxArray *outputSet[],
				 int inputs, const mxArray *inputSet[])
{
	enum
	{
		xIndex,
		maxRelativeErrorIndex,
		kNearestIndex,
		threadsIndex
	};

	std::vector<SignalPtr> xEnsemble;
	getSignals(inputSet[xIndex], std::back_inserter(xEnsemble));

	const real maxRelativeError = getReal(inputSet[maxRelativeErrorIndex]);
	const integer kNearest = getInteger(inputSet[kNearestIndex]);
	const integer threads = getInteger(inputSet[threadsIndex]);
	setNumberOfThreads(threads);
	
	outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	*rawResult = differentialEntropyKl(
		randomAccessRange(xEnsemble.begin(), xEnsemble.end()), 
		maxRelativeError, kNearest);
}
