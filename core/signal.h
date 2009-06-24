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

	class TIMCORE Signal;
	typedef Pastel::CountedPtr<Signal> SignalPtr;

	//! A finite sequence in R^n.
	/*!
	A signal is a finite sequence (m-tuple) in R^n. 
	Its dimension is n. Its size is m. The data
	of a signal can be accessed in two ways:
	via a 2-dimensional nxm array view, or
	as a sequence of points in R^n.
	Both refer to the same memory and thus can
	be used to modify the signal.
	*/

	class TIMCORE Signal
		: public Pastel::ReferenceCounted
	{
	public:
		//! Constructs a zero signal.
		/*!
		Preconditions:
		dimension >= 0
		size >= 0

		Time complexity:
		O(dimension * size)
		*/
		Signal(
			integer dimension,
			integer size,
			real defaultData = 0);

		//! Constructs an aliasing signal.
		/*!
		Preconditions:
		dimension >= 0
		size >= 0
		dimensionBegin >= 0
		dimensionBegin < signalToAlias->dimension()
		dimensionEnd >= 0
		dimensionEnd <= signalToAlias->dimension()

		Time complexity:
		constant

		An aliasing signal is one which shares
		memory with another signal.
		*/
		Signal(
			const SignalPtr& signalToAlias,
			integer dimension,
			integer dimensionBegin =  0);

		//! Swaps contents with another signal.
		/*!
		Time complexity:
		constant
		*/
		void swap(Signal& that);

		//! Returns the dimension of the signal.
		/*!
		Time complexity:
		constant
		*/
		integer dimension() const;

		//! Returns the length of the signal.
		/*!
		Time complexity:
		constant
		*/
		integer size() const;

		//! Returns a 2d view to the data.
		/*!
		Time complexity:
		constant

		The x-axis is for elements, the y-axis
		is for dimensions.
		*/
		SignalView view();

		//! Returns a non-mutable 2d view to the data.
		/*!
		Time complexity:
		constant

		The x-axis is for elements, the y-axis
		is for dimensions.
		*/
		ConstSignalView constView() const;

		//! Returns the elements as multi-dimensional points.
		/*!
		Time complexity:
		constant
		*/
		const std::vector<DynamicPoint>& pointSet() const;

		//! Returns the index:th element of the signal.
		/*!
		Preconditions:
		index >= 0
		index < size()

		Time complexity:
		constant
		*/
		DynamicPoint& operator[](integer index);

		//! Returns the index:th element of the signal.
		/*!
		Preconditions:
		index >= 0
		index < size()

		Time complexity:
		constant
		*/
		const DynamicPoint& operator[](integer index) const;

	private:
		void updatePointSet();

		Array<2, real> data_;
		std::vector<DynamicPoint> pointSet_;
		integer dimensionBegin_;
		integer dimension_;
	};

}

#include "tim/core/signal.hpp"

#endif
