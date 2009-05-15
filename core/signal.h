#ifndef TIMCORE_SIGNAL_H
#define TIMCORE_SIGNAL_H

#include "tim/core/mytypes.h"

#include <pastel/sys/countedptr.h>
#include <pastel/sys/array.h>
#include <pastel/sys/arrayview.h>

#include <vector>

namespace Tim
{

	class TIMCORE Signal
		: public Pastel::ReferenceCounted
	{
	public:
		Signal(
			integer dimension,
			integer size);

		integer dimension() const;
		integer size() const;

		View<2, real, ArrayView<2, Array<2, real> > > view();
		ConstView<2, real, ConstArrayView<2, Array<2, real> > > constView() const;

	private:
		integer dimension_;
		Array<2, real> data_;
	};

	typedef Pastel::CountedPtr<Signal> SignalPtr;

	SignalPtr mergeSignalDimensions(
		const std::vector<SignalPtr>& signalList);

}

#endif
