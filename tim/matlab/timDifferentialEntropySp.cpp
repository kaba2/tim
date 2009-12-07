#include "tim/matlab/tim_mex.h"

#include "tim/core/differential_entropy_sp.h"

using namespace Tim;

void mexFunction(int outputs, mxArray *outputSet[],
				 int inputs, const mxArray *inputSet[])
{
	enum
	{
		xIndex
	};

	std::vector<SignalPtr> xEnsemble;
	getSignals(inputSet[xIndex], std::back_inserter(xEnsemble));

	outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	*rawResult = differentialEntropySp(
		forwardRange(xEnsemble.begin(), xEnsemble.end()));
}
