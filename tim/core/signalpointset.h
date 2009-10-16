// Description: SignalPointSet class
// Detail: Models a semi-dynamic point set

#ifndef TIM_SIGNALPOINTSET_H
#define TIM_SIGNALPOINTSET_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include <pastel/geometry/pointkdtree.h>

#include <pastel/sys/countedptr.h>
#include <pastel/sys/forwardrange.h>
#include <pastel/sys/array.h>

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
		// The Array_ObjectPolicy_PointKdTree<real> represents points
		// by a const real* to the beginning of point coordinate data.
		typedef PointKdTree<real, Dynamic, Array_ObjectPolicy_PointKdTree<real> > KdTree;
		typedef KdTree::ConstObjectIterator ConstObjectIterator;
		typedef KdTree::Object Object;

	private:
		// This container will hold the set of objects currently in
		// the kd-tree and thus also in the time-window. Because
		// objects need to be inserted and removed from both front
		// and back, the std::deque is a good choice for a container.
		typedef std::deque<ConstObjectIterator> ObjectSet;

	public:
		// Using default destructor.

		typedef ObjectSet::const_iterator ConstObjectIterator_Iterator;

		//! Constructs using the given ensemble of signals.
		template <typename SignalPtr_Iterator>
		explicit SignalPointSet(
			const ForwardRange<SignalPtr_Iterator>& signalSet,
			bool startFull = false);

		//! Constructs using given subdimensions and initial time-window.
		/*!
		Preconditions:
		dimensionBegin <= dimensionEnd
		dimensionBegin >= 0
		dimensionEnd <= signalSet.front()->dimension()

		The subdimension integer interval is
		given by [dimensionBegin, dimensionEnd[. If 'startFull' is true,
		then initially all samples are contained in the time-window.
		Otherwise no samples are initially contained in the time-window.
		*/
		template <typename SignalPtr_Iterator>
		SignalPointSet(
			const ForwardRange<SignalPtr_Iterator>& signalSet,
			bool startFull,
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

		//! First iterator to set of points currently in the kd-tree.
		ConstObjectIterator_Iterator begin() const;

		//! One-past-last iterator to set of points currently in the kd-tree.
		ConstObjectIterator_Iterator end() const;

		//! Returns the beginning time of the current time-window.
		/*!
		The interval of the time-window is given by	[timeBegin(), timeEnd()[.
		*/
		integer timeBegin() const;

		//! Returns the one-past-last time of the current time-window.
		/*!
		The interval of the time-window is given by	[timeBegin(), timeEnd()[.
		*/
		integer timeEnd() const;

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

		//! Returns a point that corresponds to a given kd-tree object.
		VectorD point(const Object& object) const;

	private:
		// Prohibited, for now.
		SignalPointSet(const SignalPointSet& that);

		// Extracts points from the given ensemble of signals.
		/*
		Each point is a pointer to the beginning of its coordinate
		data. These pointers are placed in 'pointSet_'.
		*/
		template <typename SignalPtr_Iterator>
		void createPointSet(
			const ForwardRange<SignalPtr_Iterator>& signalSet);

		void construct(bool startFull);

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
		time-window. Note that kdTree_ uses an object policy
		in which the object is a const real*. This container is
		needed to insert points when the time-window is moved.

		objectSet_:
		Contains a set of object iterators to the kdTree_,
		including only those points which are currently in
		the time-window. This container is needed to be able
		to remove points when the time-window is moved.

		samples_:
		Contains the number of samples that are considered
		for each signal.

		timeBegin_, timeEnd_:
		The position of the time-window is given by the
		integer interval [timeBegin_, timeEnd_[.

		dimensionBegin_, dimension_:
		The signals can also be considered by their
		marginal signals. The integer interval 
		[dimensionBegin_, dimensionBegin_ + dimension_]
		denotes a subdimension interval inside
		the signals.
		*/

		KdTree kdTree_;
		std::vector<const real*> pointSet_;
		ObjectSet objectSet_;
		integer signals_;
		integer samples_;
		integer timeBegin_;
		integer timeEnd_;
		integer dimensionBegin_;
		integer dimension_;
	};

	typedef CountedPtr<SignalPointSet> SignalPointSetPtr;

}

#include "tim/core/signalpointset.hpp"

#endif
