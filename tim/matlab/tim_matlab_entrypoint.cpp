#include "tim/matlab/tim_matlab.h"

// This is the actual mex entry-point function.

void mexFunction(
	int outputs, mxArray *outputSet[],
	int inputs, const mxArray *inputSet[])
{
	Tim::matlabCallFunction(
		outputs, outputSet,
		inputs, inputSet);
}
