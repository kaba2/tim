// Documentation: tim_matlab_cpp.txt

#include "tim/matlab/tim_matlab.h"

#include <map>

namespace Tim
{

	typedef std::map<std::string, MatlabFunction*> FunctionMap;
	typedef FunctionMap::const_iterator FunctionIterator;	

	FunctionMap functionMap;

	void matlabAddFunction(
		const std::string& name,
		MatlabFunction* function)
	{
		functionMap.insert(std::make_pair(name, function));
	}

}

using namespace Tim;

namespace
{

	void matlabListFunctions(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		// This function lists all the functions
		// that are callable with tim_matlab().
		// It is itself such a function.

		FunctionIterator iter = functionMap.begin();
		const FunctionIterator iterEnd = functionMap.end();
		
		while(iter != iterEnd)
		{
			printf(iter->first.c_str());
			printf("\n");
			++iter;
		}
	}

	void addFunction()
	{
		matlabAddFunction(
			"list_functions",
			matlabListFunctions);
	}

	CallFunction run(addFunction);

}

// Actual Mex entry-point for Matlab.

void mexFunction(
	int outputs, mxArray *outputSet[],
	int inputs, const mxArray *inputSet[])
{
	// The first parameter is the name of the
	// function inside TIM Matlab that should
	// be called.

	const std::string name = getString(inputSet[0]);

	// See if a function with that name has
	// been registered.
	
	FunctionIterator iter = functionMap.find(name);
	
	if (iter != functionMap.end())
	{
		// There is a function with that
		// name registered.

		// Find out its function pointer.
		MatlabFunction* function = iter->second;

		// Call that function, but omit the
		// first name parameter.
		(*function)(
			outputs, outputSet,
			inputs - 1, inputSet + 1);
	}
}