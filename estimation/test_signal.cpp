#include "estimation.h"

#include "tim/core/signal_tools.h"

#include <pastel/sys/view_all.h>

using namespace Tim;

namespace
{

	void testSignal()
	{
		// Create a signal of dimension 4
		// and of size 10.

		const SignalPtr signal = 
			SignalPtr(new Signal(10, 4));

		// The signals can be output to a stream
		// to see the contents.

		std::cout << *signal << std::endl;

		const integer dimension = signal->width();
		const integer samples = signal->height();

		// The most primitive way to view the signal data
		// is as a collection of numbers. One obtains
		// the element at (y, x)

		// One way to view the signal data is as
		// a collection of samples. One obtains
		// the y:th point by the bracket notation [y].
		// See 'pastel/sys/point.h' for info
		// on samples.

		for (integer y = 0;y < samples;++y)
		{
			for (integer x = 0;x < dimension;++x)
			{
				(*signal)[y][x] = x * y;
			}
		}

		std::cout << *signal << std::endl;

		// Another way to view the signal data is
		// as an 2d array of values. Such arrays
		// can be handled efficiently via so-called
		// views. A view is an abstract array
		// (in generic programming sense), which
		// allows to concentrate on subsets of
		// a concrete array without copying it.
		// Because of being based on generic programming
		// and cursors (similar to iterators), 
		// the resulting code should be close to
		// performance to an equivalent hand-written code.
		// For example, one could take
		// some rectangular subset of the array, and
		// then out of that subset pick only those
		// coordinates which are odd, and finally
		// assign -1 to those elements.
		// Let us do that.

		clear(-1, 
			sparseView(
			subView(signal->view(), Rectangle2(0, 3, 4, 7)), 
			Point2i(1, 1), Vector2i(2, 2)));

		std::cout << *signal << std::endl;
	}

	void addTest()
	{
		timTestList().add("signal", testSignal);
	}

	CallFunction run(addTest);

}

