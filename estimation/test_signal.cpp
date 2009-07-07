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

		// The signals can be output to a stream
		// to see the contents.

		std::cout << *signal << std::endl;

		const integer dimension = signal->dimension();
		const integer samples = signal->samples();

		// The most primitive way to view the signal data
		// is as a collection of numbers. One obtains
		// the element at (y, x)

		// One way to view the signal data is as
		// a collection of samples. One obtains
		// the y:th point by the bracket notation [y].
		// See 'pastel/sys/point.h' for info
		// on samples.

		for (integer x = 0;x < dimension;++x)
		{
			for (integer y = 0;y < samples;++y)
			{
				signal->data()[x][y] = x * y;
			}
		}

		std::cout << *signal << std::endl;
	}

	void addTest()
	{
		timTestList().add("signal", testSignal);
	}

	CallFunction run(addTest);

}

