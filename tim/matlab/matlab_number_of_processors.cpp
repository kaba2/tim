// Description: number_of_processors
// DocumentationOf: number_of_processors.m

#include "tim/matlab/tim_matlab.h"

void force_linking_number_of_processors() {};

using namespace Tim;

namespace
{

	void matlabNumberOfProcessors(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum
		{
			Inputs
		};

		enum Output
		{
			Processors,
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		integer* processors = createScalar<integer>(
			outputSet[Processors]);
		
		*processors = numberOfProcessors();
	}

	void addFunction()
	{
		matlabAddFunction(
			"number_of_processors",
			matlabNumberOfProcessors);
	}

	CallFunction run(addFunction);

}
