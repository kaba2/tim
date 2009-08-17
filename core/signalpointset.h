#ifndef TIM_SIGNALPOINTSET_H
#define TIM_SIGNALPOINTSET_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include <pastel/geometry/pointkdtree.h>

#include <pastel/sys/countedptr.h>
#include <pastel/sys/forwardrange.h>
#include <pastel/sys/array.h>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/reverse_iterator.hpp>

#include <deque>

namespace Tim
{

	class TIMCORE SignalPointSet_TimeWindow
	{
	public:
		enum Enum
		{
			StartFull,
			StartEmpty
		};
	};

	class TIMCORE SignalPointSet
		: public ReferenceCounted
	{
	public:
		typedef PointKdTree<real, Dynamic, Array_ObjectPolicy_PointKdTree<real> > KdTree;
		typedef KdTree::ConstObjectIterator ConstObjectIterator;

	private:
		typedef std::deque<ConstObjectIterator> ObjectSet;

	public:
		// Using default destructor.

		typedef ObjectSet::const_iterator ConstObjectIterator_Iterator;

		template <typename Signal_Iterator>
		explicit SignalPointSet(
			const ForwardRange<Signal_Iterator>& signalSet,
			SignalPointSet_TimeWindow::Enum timeWindowStart = 
			SignalPointSet_TimeWindow::StartEmpty);

		template <typename Signal_Iterator>
		explicit SignalPointSet(
			const ForwardRange<Signal_Iterator>& signalSet,
			SignalPointSet_TimeWindow::Enum timeWindowStart,
			integer dimensionBegin,
			integer dimensionEnd);

		//! Swaps the contents of two SignalPointSet's.
		void swap(SignalPointSet& that);

		//! Copies the contents of another SignalPointSet.
		SignalPointSet& operator=(const SignalPointSet& that);

		//! Set the time window.
		/*!
		Preconditions:
		tNewBegin <= tNewEnd
		tNewBegin >= 0
		tNewEnd <= samples()
		*/
		void setTimeWindow(integer tNewBegin, integer tNewEnd);

		const KdTree& kdTree() const;

		ConstObjectIterator_Iterator begin() const;
		ConstObjectIterator_Iterator end() const;

		integer timeBegin() const;
		integer timeEnd() const;

		integer samples() const;

		integer dimension() const;

	private:
		// Prohibited, for now.
		SignalPointSet(const SignalPointSet& that);

		template <typename Signal_Iterator>
		void construct(
			const ForwardRange<Signal_Iterator>& signalSet,
			SignalPointSet_TimeWindow::Enum timeWindowStart,
			integer dimensionBegin,
			integer dimensionEnd);

		KdTree kdTree_;
		std::vector<SignalPtr> signalSet_;
		ObjectSet objectSet_;
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
