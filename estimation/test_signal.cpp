#include "estimation.h"

#include "tim/core/signal_tools.h"

#include <pastel/sys/view_all.h>
#include <pastel/sys/random.h>
#include <pastel/sys/string_tools.h>

#include "tim/core/mutual_information.h"

using namespace Tim;

namespace
{

	void testSignal()
	{
		// Create a signal of dimension 4
		// and of size 10.

		const SignalPtr signal = 
			SignalPtr(new Signal(4, 10));

		signal->setName("Experiment 1");

		// The signals can be output to a stream
		// to see the contents.

		std::cout << *signal << std::endl;

		const integer dimension = signal->dimension();
		const integer samples = signal->samples();

		// The most primitive way to view the signal data
		// is as a collection of numbers. One obtains
		// the element at (y, x)

		for (integer i = 0;i < dimension;++i)
		{
			for (integer j = 0;j < samples;++j)
			{
				signal->data()(i, j) = i * j;
			}
		}

		std::cout << *signal << std::endl;
	}

	void testSplit()
	{
		SignalPtr xy = 
			SignalPtr(new Signal(2, 5));

		SignalPtr x = slice(xy, 0, 2);
		SignalPtr y = slice(xy, 2, 5);
	}

	void addTest()
	{
		timTestList().add("signal", testSignal);
	}

	CallFunction run(addTest);

}

