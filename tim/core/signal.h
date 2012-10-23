// Description: Signal class
// Detail: Models a time series

#ifndef TIM_SIGNAL_H
#define TIM_SIGNAL_H

#include <pastel/math/matrix.h>

#include <pastel/sys/countedptr.h>
#include <pastel/sys/sparse_iterator.h>
#include <pastel/sys/counting_iterator.h>

#include "tim/core/mytypes.h"

#include <vector>

namespace Tim
{

	//! A time series.

	class TIM Signal
		: public ReferenceCounted
	{
	private:
		typedef SparseIterator<CountingIterator<real*> > 
			Point_Iterator;
		typedef ConstSparseIterator<CountingIterator<const real*> > 
			Point_ConstIterator;

	public:
		// Using default copy constructor.
		// Using default assignment.
		// Using default destructor.

		//! Construct a signal with 0 samples and 0 dimension.
		Signal();

		//! Construct a zero signal.
		/*!
		Preconditions:
		samples >= 0
		dimension >= 0
		*/
		Signal(integer samples, integer dimension, integer t = 0);

		//! Construct an aliasing signal.
		/*!
		Preconditions:
		samples >= 0
		dimension >= 0

		Using aliasing, the signal will refer to the given
		memory region, rather than allocating its own.
		No attempts will be made to deallocate that memory.
		This can be used to efficiently form sub-signals
		which refer to a sub-ranges of dimensions for a
		given signal.
		*/
		Signal(integer samples, integer dimension, 
			integer t, real* dataToAlias);

		//! Sets a name for the signal.
		void setName(const std::string& name);

		//! Returns the name of the signal.
		const std::string& name() const;

		//! Returns the dimension of the signal.
		integer dimension() const;

		//! Returns the number of samples in the signal.
		integer samples() const;

		//! Returns the samples in a matrix.
		/*!
		Each row contains one sample.
		*/
		Matrix<real>& data();

		//! Returns the samples in a matrix.
		/*!
		See the documentation for the non-const version.
		*/
		const Matrix<real>& data() const;

		//! Sets the time position of the first sample.
		void setT(integer t);

		//! Returns the time position of the first sample.
		integer t() const;

		//! Returns an iterator to the beginning of the point set.
		/*!
		Here the signal data is interpreted so that each sample 
		is a point which can be accessed by the Array_PointPolicy<real>
		point policy. See 'pastel/sys/pointpolicy.txt'.
		*/
		Point_Iterator pointBegin(integer dimensionBegin = 0);

		//! Returns an iterator to the beginning of the point set.
		/*!
		See the documentation for the non-const version.
		*/
		Point_ConstIterator pointBegin(integer dimensionBegin = 0) const;

		//! Returns an iterator to the end of the point set.
		/*!
		See the documentation for 'pointBegin'.
		*/
		Point_Iterator pointEnd(integer dimensionBegin = 0);

		//! Returns an iterator to the end of the point set.
		/*!
		See the documentation for 'pointBegin'.
		*/
		Point_ConstIterator pointEnd(integer dimensionBegin = 0) const;

	private:
		//! The name of the signal.
		std::string name_;
		
		//! The signal data.
		/*!
		Each row is one sample, and thus the width
		of this matrix is the dimensionality. The number
		of samples is given by the number of rows.
		*/
		Matrix<real> data_;
		
		//! The time position of the first sample.
		integer t_;
	};

	typedef CountedPtr<Signal> SignalPtr;

}

#include "tim/core/signal.hpp"

#endif
