#ifndef TIM_SIGNAL_H
#define TIM_SIGNAL_H

#include "tim/core/mytypes.h"

#include <pastel/sys/countedptr.h>
#include <pastel/sys/array.h>
#include <pastel/sys/arrayview.h>

#include <vector>

namespace Tim
{

	typedef View<2, real, ArrayView<2, Array<2, real> > > SignalView;
	typedef ConstView<2, real, ConstArrayView<2, Array<2, real> > > ConstSignalView;

	class TIMCORE Signal
		: public Pastel::ReferenceCounted
	{
	public:
		Signal(
			integer dimension,
			integer size);

		integer dimension() const;
		integer size() const;

		SignalView view();
		ConstSignalView constView() const;

		const std::vector<DynamicPoint>& pointSet() const;

	private:
		integer dimension_;
		Array<2, real> data_;
		std::vector<DynamicPoint> pointSet_;
	};

	typedef Pastel::CountedPtr<Signal> SignalPtr;

	TIMCORE SignalPtr mergeSignalDimensions(
		const std::vector<SignalPtr>& signalList);

	SignalPtr newSignal(integer dimension, integer size);

}

#include "tim/core/signal.hpp"

#endif
