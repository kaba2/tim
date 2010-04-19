// Description: Signal class
// Detail: Models a time series

#ifndef TIM_SIGNAL_H
#define TIM_SIGNAL_H

#include "tim/core/mytypes.h"

#include <pastel/math/matrix.h>

#include <pastel/sys/countedptr.h>
#include <pastel/sys/sparseiterator.h>
#include <pastel/sys/countingiterator.h>

#include <vector>

namespace Tim
{

	//! A time series.

	class TIM Signal
		: public ReferenceCounted
	{
	private:
		typedef SparseIterator<CountingIterator<real*> > PointIterator;
		typedef ConstSparseIterator<CountingIterator<const real*> > ConstPointIterator;

	public:
		// Using default copy constructor.
		// Using default assignment.
		// Using default destructor.

		Signal();
		Signal(integer samples, integer dimension, integer t = 0);
		Signal(integer samples, integer dimension, 
			integer t, real* dataToAlias);

		void setName(const std::string& name);
		const std::string& name() const;

		integer dimension() const;
		integer samples() const;

		MatrixD& data();
		const MatrixD& data() const;

		void setT(integer t);
		integer t() const;

		PointIterator pointBegin(integer dimensionBegin = 0);
		ConstPointIterator pointBegin(integer dimensionBegin = 0) const;

		PointIterator pointEnd(integer dimensionBegin = 0);
		ConstPointIterator pointEnd(integer dimensionBegin = 0) const;

	private:
		std::string name_;
		MatrixD data_;
		integer t_;
	};

	typedef CountedPtr<Signal> SignalPtr;

}

#include "tim/core/signal.hpp"

#endif
