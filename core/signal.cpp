#include "tim/core/signal.h"

#include <pastel/sys/subview.h>
#include <pastel/sys/view_tools.h>

namespace Tim
{

	Signal::Signal(integer dimension, integer size)
		: dimension_(dimension)
		, data_(dimension, size, 0)
		, pointSet_()
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(size >= 0, size);

		pointSet_.resize(size, nullPoint<Dynamic, real>());

		for (integer i = 0;i < size;++i)
		{
			DynamicPoint temp(
				aliasPoint<Dynamic, real>(dimension, &data_(0, i)));

			pointSet_[i].swap(temp);

			ASSERT(&pointSet_[i][0] == &data_(0, i));
		}
	}

	integer Signal::size() const
	{
		return data_.height();
	}

	integer Signal::dimension() const
	{
		return data_.width();
	}

	SignalView Signal::view()
	{
		return arrayView(data_);
	}

	ConstSignalView Signal::constView() const
	{
		return constArrayView(data_);
	}

	const std::vector<DynamicPoint>& Signal::pointSet() const
	{
		return pointSet_;
	}

	TIMCORE SignalPtr mergeSignalDimensions(
		const std::vector<SignalPtr>& signalList)
	{
		if (signalList.empty())
		{
			return SignalPtr();
		}

		const integer size = signalList.front()->size();

		const integer signals = signalList.size();
		integer bigDimension = 0;
		for (integer i = 0;i < signals;++i)
		{
			bigDimension += signalList[i]->dimension();
			ENSURE2(signalList[0]->size() == signalList[i]->size(), 
				signalList[0]->size(), signalList[i]->size());
		}

		SignalPtr bigSignal(new Signal(bigDimension, size));
		
		integer dimensionOffset = 0;

		for (integer i = 0;i < signals;++i)
		{
			copy(signalList[i]->constView(),
				subView(bigSignal->view(), 
				Rectangle2(dimensionOffset, 0, 
				dimensionOffset + signalList[i]->dimension(), size)));

			dimensionOffset += signalList[i]->dimension();
		}

		return bigSignal;
	}

}
