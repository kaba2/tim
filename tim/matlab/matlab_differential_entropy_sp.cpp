// Description: differential_entropy_sp
// Documentation: tim_matlab_functions.txt

#include "tim/matlab/tim_matlab.h"

#include "tim/core/differential_entropy_sp.h"

using namespace Tim;

namespace
{

	void matlabDifferentialEntropySp(
		int outputs, mxArray *outputSet[],
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
			range(xEnsemble.begin(), xEnsemble.end()));
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_sp",
			matlabDifferentialEntropySp);
	}

	CallFunction run(addFunction);

}
