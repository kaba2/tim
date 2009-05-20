#include "estimation.h"

namespace Tim
{

	void testMutualInformation()
	{
		
	}

	void testAdd()
	{
		timTestList().add("MutualInformation", testMutualInformation);
	}

	CallFunction run(testAdd);

}