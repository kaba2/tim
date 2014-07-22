// Description: SignalPointSet class
// Detail: Models a semi-dynamic point set

#ifndef TIM_SIGNALPOINTSET_H
#define TIM_SIGNALPOINTSET_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include <pastel/geometry/pointkdtree.h>

#include <pastel/sys/countedptr.h>
#include <pastel/sys/range.h>
#include <pastel/sys/array.h>
#include <pastel/sys/pointer_locator.h>

#include <vector>
#include <deque>

namespace Tim
{

	//! SignalPointSet class
	/*!
	Models a semi-dynamic point set. It is reference counted
	using the CountedPtr smart pointer from the Pastel library.
	*/
	class TIM SignalPointSet
		: public ReferenceCounted
	{
	public:
		// The Array_PointPolicy<real> represents points
		// by a const real* to the beginning of point coordinate data.
		using Settings = PointKdTree_Settings<Pointer_Locator<real>>;
		typedef PointKdTree<Settings> KdTree;
		typedef KdTree::Point_ConstIterator Point_ConstIterator;
		typedef KdTree::Point Point;

		typedef std::vector<Point_ConstIterator> PointSet;
		typedef PointSet::const_iterator Point_ConstIterator_Iterator;

	public:
		// Using default destructor.

		//! Constructs using the given ensemble of signals.
		template <typename Signal_Iterator>
		explicit SignalPointSet(
			const boost::iterator_range<Signal_Iterator>& signalSet);

		//! Constructs using given subdimensions and initial time-window.
		/*!
		Preconditions:
		dimensionBegin <= dimensionEnd
		dimensionBegin >= 0
		dimensionEnd <= signalSet.front().dimension()

		The subdimension integer interval is
		given by [dimensionBegin, dimensionEnd[. If 'startFull' is true,
		then initially all samples are contained in the time-window.
		Otherwise no samples are initially contained in the time-window.
		*/
		template <typename Signal_Iterator>
		SignalPointSet(
			const boost::iterator_range<Signal_Iterator>& signalSet,
			integer dimensionBegin,
			integer dimensionEnd);

		//! Swaps the contents of two SignalPointSet's.
		/*!
		Time complexity: constant
		*/
		void swap(SignalPointSet& that);

		//! Copies the contents of another SignalPointSet.
		SignalPointSet& operator=(const SignalPointSet& that);

		//! Set the time-window.
		/*!
		Preconditions:
		tNewBegin <= tNewEnd

		The more the old time-window and the new time-window
		overlap, the less work needs to be done.
		*/
		void setTimeWindow(integer tNewBegin, integer tNewEnd);

		//! Returns a non-mutable reference to the multi-resolution kd-tree.
		const KdTree& kdTree() const;

		//! First iterator to set of points currently in the window.
		Point_ConstIterator_Iterator begin() const;

		//! One-past-last iterator to set of points currently in the window.
		Point_ConstIterator_Iterator end() const;

		//! Returns the beginning time of the current time-window.
		/*!
		The interval of the time-window is given by	[windowBegin(), windowEnd()[.
		*/
		integer windowBegin() const;

		//! Returns the one-past-last time of the current time-window.
		/*!
		The interval of the time-window is given by	[windowBegin(), windowEnd()[.
		*/
		integer windowEnd() const;

		//! Returns the total number of points in the point set.
		integer samples() const;

		//! Returns the dimension of the point set.
		/*!
		Note that the dimension of the point set is not
		necessarily the same as the dimension of the input
		signals. This is because you can choose a subinterval
		of the dimensions in construction.
		*/
		integer dimension() const;

		//! Returns the dimension offset of the signal data.
		/*!
		Using a dimension offset allows you to select a subrange
		of the dimensions of a signal.
		*/
		integer dimensionBegin() const;

		//! Returns a vector that corresponds to a given point.
		VectorD point(const Point& object) const;

	private:
		// Prohibited, for now.
		SignalPointSet(const SignalPointSet& that);

		// Extracts points from the given ensemble of signals.
		/*
		Each point is a pointer to the beginning of its coordinate
		data. These pointers are placed in 'pointSet_'.
		*/
		template <typename Signal_Iterator>
		void createPointSet(
			const boost::iterator_range<Signal_Iterator>& signalSet);

		void hide(
			const AlignedBox<integer, 1>& range);

		void show(
			const AlignedBox<integer, 1>& range);

		/*
		kdTree_:
		A multi-resolution kd-tree that has been subdivided with
		all available points but which only contains the points
		in the time-window at any time instant. Note multi-resolution
		kd-tree adapts itself automatically to smaller number of points.

		signalSet_:
		Contains the ensemble of signals that are the input 
		data.

		pointSet_:
		Contains a set of pointers with each pointer pointing
		to the beginning of point coordinate data. This set represents
		the set of all available points without culling by a 
		time-window. Note that kdTree_ uses a PointPolicy
		in which the point is a const real*. This container is
		needed to insert points when the time-window is moved.

		samples_:
		Contains the number of samples that are considered
		for each signal.

		windowBegin_, windowEnd_:
		The position of the time-window is given by the
		integer interval [windowBegin_, windowEnd_[.

		dimensionBegin_, dimension_:
		The signals can also be considered by their
		marginal signals. The integer interval 
		[dimensionBegin_, dimensionBegin_ + dimension_]
		denotes a subdimension interval inside
		the signals.

		timeBegin_:
		The time instant t corresponding to 'pointSet_[i]'
		is given by 't = timeBegin_ + (i / signals_)'.
		*/

		KdTree kdTree_;
		PointSet pointSet_;
		integer signals_;
		integer samples_;
		integer windowBegin_;
		integer windowEnd_;
		integer dimensionBegin_;
		integer dimension_;
		integer timeBegin_;
	};

	typedef CountedPtr<SignalPointSet> SignalPointSetPtr;

}

#include "tim/core/signalpointset.hpp"

#endif
