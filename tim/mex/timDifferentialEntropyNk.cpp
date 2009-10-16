#include "tim/mex/tim_mex.h"

#include "tim/core/differential_entropy_nk.h"

#include <pastel/math/euclidean_normbijection.h>

using namespace Tim;

void mexFunction(int outputs, mxArray *outputSet[],
				 int inputs, const mxArray *inputSet[])
{
	enum
	{
		xIndex,
		maxRelativeErrorIndex,
		threadsIndex
	};

	std::vector<SignalPtr> xEnsemble;
	getSignals(inputSet[xIndex], std::back_inserter(xEnsemble));

	const real maxRelativeError = getReal(inputSet[maxRelativeErrorIndex]);
	const integer threads = getInteger(inputSet[threadsIndex]);
	setNumberOfThreads(threads);

	integer intrinsicDimension = 0;

	const real entropy = 
		differentialEntropyNk(
		forwardRange(xEnsemble.begin(), xEnsemble.end()), 
		maxRelativeError,
		Euclidean_NormBijection<real>(),
		&intrinsicDimension);

	if (outputs > 0)
	{
		outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
		real* outEntropy = mxGetPr(outputSet[0]);
		*outEntropy = entropy;
	}

	if (outputs > 1)
	{
		outputSet[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
		real* outIntrinsicDimension = mxGetPr(outputSet[1]);
		*outIntrinsicDimension = intrinsicDimension;
	}
}
