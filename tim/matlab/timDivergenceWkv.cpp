#include "tim/matlab/tim_mex.h"

#include "tim/core/divergence_wkv.h"

using namespace Tim;

void mexFunction(int outputs, mxArray *outputSet[],
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

	const integer threads = getInteger(inputSet[threadsIndex]);
	setNumberOfThreads(threads);

	outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	*rawResult = divergenceWkv(
		randomAccessRange(xEnsemble.begin(), xEnsemble.end()), 
		randomAccessRange(yEnsemble.begin(), yEnsemble.end()));
}
