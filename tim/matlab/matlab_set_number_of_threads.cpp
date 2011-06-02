// Description: set_number_of_threads
// DocumentationOf: set_number_of_threads.m

#include "tim/matlab/tim_matlab.h"

void force_linking_set_number_of_threads() {};

using namespace Tim;

namespace
{

	void matlabSetNumberOfThreads(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum
		{
			Threads,
			Inputs
		};

		enum Output
		{
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		const integer threads = 
			asScalar<integer>(inputSet[Threads]);

		ENSURE_OP(threads, >, 0);
		
		setNumberOfThreads(threads);
	}

	void addFunction()
	{
		matlabAddFunction(
			"set_number_of_threads",
			matlabSetNumberOfThreads);
	}

	CallFunction run(addFunction);

}
