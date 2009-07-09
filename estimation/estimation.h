#ifndef TIM_ESTIMATION_H
#define TIM_ESTIMATION_H

#include <pastel/sys/testrunner.h>
#include <pastel/sys/callfunction.h>
#include <pastel/sys/testsuite.h>

#include "tim/core/mytypes.h"

namespace Tim
{

	inline Pastel::TestRunner& timTestList()
	{
		static Pastel::TestRunner timTestRunner("Tim library");
		return timTestRunner;
	}

	inline Pastel::TestReport& timTestReport()
	{
		static Pastel::TestReport theTestReport("Tim");
		return theTestReport;
	}

}

#endif
