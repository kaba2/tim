#ifndef TIM_ESTIMATION_H
#define TIM_ESTIMATION_H

#include <pastel/sys/testrunner.h>
#include <pastel/sys/callfunction.h>

#include "tim/core/mytypes.h"

namespace Tim
{

	inline Pastel::TestRunner& timTestList()
	{
		static Pastel::TestRunner timTestRunner("Tim library");
		return timTestRunner;
	}


}

#endif
