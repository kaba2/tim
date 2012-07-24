#include "tim/core/signal_tools.h"

#include <pastel/sys/countingiterator.h>
#include <pastel/sys/sparseiterator.h>
#include <pastel/sys/nulliterator.h>

#include <vector>
#include <list>
#include <set>

using namespace Tim;

namespace
{

	void basicUseCase()
	{
		// Create a signal.

		SignalPtr xy = 
			SignalPtr(new Signal(5, 3));

		// The data() function returns a reference
		// to the data matrix. This allows to use
		// all the tools available for matrices. Here
		// we use the ability to assign values
		// via a comma separated list.

		xy->data() |=
			0, 5, 10,
			1, 6, 11,
			2, 7, 12,
			3, 8, 13,
			4, 9, 14;
	}

	void splitUseCase()
	{
		SignalPtr xy = SignalPtr(new Signal(100, 3));

		// Given a signal, the split() function
		// can be used to select a range of dimensions
		// which are used to form a new signal.

		SignalPtr x = split(xy, 0, 1);
		SignalPtr y = split(xy, 1, 3);

		// The following split() function
		// outputs all the 1d subsignals of the given
		// signal. Many of the functions output their 
		// results through an output iterator. This also
		// works as an example how iterators can be used 
		// to output to a container. 

		std::vector<SignalPtr> splitSet;
		split(xy, std::back_inserter(splitSet));

		// The split() function can be given a partition
		// set which determines the dimensions of the
		// subsignals to be splitted. Here the subsignals
		// will cover dimensions [0, 1[ and [1, 3[:

		std::vector<integer> partition;
		partition.push_back(0);
		partition.push_back(1);
		partition.push_back(3);
		split(xy, partition, std::back_inserter(splitSet));
	}

	void mergeUseCase()
	{
		SignalPtr x;
		SignalPtr y;

		// Signals can be merged by using the
		// merge() function:

		SignalPtr yx = merge(y, x);
		
		// Let us pack signals into a container.

		std::vector<SignalPtr> signalSet;
		signalSet.push_back(x);
		signalSet.push_back(y);
		signalSet.push_back(yx);
		
		// merge() can work with sets of signals by 
		// the iterator range abstraction:

		SignalPtr xyyx = merge(
			range(signalSet.begin(), signalSet.end()));

		// You could have also used any other container:

		std::set<SignalPtr> ySignalSet;
		ySignalSet.insert(x);
		ySignalSet.insert(y);
		ySignalSet.insert(yx);

		SignalPtr bXyyx = merge(
			range(ySignalSet.begin(), ySignalSet.end()));
		
		// When merging signals, each of the signals
		// can be delayed. In the following we delay
		// x by 3, y by 4, and yx by 5.
		
		std::vector<integer> lagSet;
		lagSet.push_back(3);
		lagSet.push_back(4);
		lagSet.push_back(5);

		SignalPtr delayedXyyx = merge(
			range(signalSet.begin(), signalSet.end()),
			range(lagSet.begin(), lagSet.end()));

		// However, since the lags are sequential number, 
		// we could have avoided forming the lagSet container
		// by using a CountingIterator instead.
		// Also note that an iterator range can be given
		// in the begin-size form.

		SignalPtr bDelayedXyyx = merge(
			range(signalSet.begin(), signalSet.end()),
			range(countingIterator(3), signalSet.size()));

		// If the lags were 3, 5, and 7, we could still use a 
		// combination of a CountingIterator and a SparseIterator:

		SignalPtr cDelayedXyyx = merge(
			range(signalSet.begin(), signalSet.end()),
			range(sparseIterator(countingIterator(3), 2), signalSet.size()));

		// Arrays can be used to form iterator ranges.

		SignalPtr signalArray[3] = {x, y, yx};
		SignalPtr eXyyx = merge(range(signalArray));
	}

}
