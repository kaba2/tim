#include "tim/core/signal.h"

#include <pastel/sys/subview.h>
#include <pastel/sys/view_tools.h>

namespace Tim
{

	Signal::Signal(integer dimension, integer size)
		: dimension_(dimension)
		, data_(dimension, size, 0)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(size >= 0, size);
	}

	integer Signal::size() const
	{
		return data_.height();
	}

	integer Signal::dimension() const
	{
		return data_.width();
	}

	View<2, real, ArrayView<2, Array<2, real> > > Signal::view()
	{
		return arrayView(data_);
	}

	ConstView<2, real, ConstArrayView<2, Array<2, real> > > Signal::constView() const
	{
		return constArrayView(data_);
	}

	SignalPtr mergeSignalDimensions(
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
