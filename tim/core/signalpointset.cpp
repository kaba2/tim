#include "tim/core/signalpointset.h"
#include "tim/core/signal_tools.h"

#include <pastel/geometry/search_all_neighbors_pointkdtree.h>
#include <pastel/geometry/slidingmidpoint_splitrule_pointkdtree.h>

#include <pastel/sys/constantiterator.h>

#include <pastel/math/normbijection.h>

#include <boost/iterator/reverse_iterator.hpp>

namespace Tim
{

	SignalPointSet::SignalPointSet(const SignalPointSet& that)
	{
		ENSURE(false);
	}

	void SignalPointSet::construct(bool startFull)
	{
		// Insert all the points into the tree.

		if (startFull)
		{
			kdTree_.insert(pointSet_.begin(), pointSet_.end(), 
				std::back_inserter(objectSet_));
			
			timeBegin_ = 0;
			timeEnd_ = pointSet_.size();
		}
		else
		{
			kdTree_.insert(pointSet_.begin(), pointSet_.end());
		}

		// Compute a fine subdivision for the points.

		kdTree_.refine(SlidingMidpoint2_SplitRule_PointKdTree());

		if (!startFull)
		{
			kdTree_.eraseObjects();
		}
	}

	void SignalPointSet::swap(SignalPointSet& that)
	{
		kdTree_.swap(that.kdTree_);
		pointSet_.swap(that.pointSet_);
		objectSet_.swap(that.objectSet_);
		std::swap(signals_, that.signals_);
		std::swap(samples_, that.samples_);
		std::swap(timeBegin_, that.timeBegin_);
		std::swap(timeEnd_, that.timeEnd_);
		std::swap(dimensionBegin_, that.dimensionBegin_);
		std::swap(dimension_, that.dimension_);
	}

	SignalPointSet& SignalPointSet::operator=(
		const SignalPointSet& that)
	{
		SignalPointSet copy(that);
		swap(copy);
		return *this;
	}

	void SignalPointSet::setTimeWindow(
		integer newTimeBegin, integer newTimeEnd)
	{
		ENSURE_OP(newTimeBegin, <=, newTimeEnd);

		newTimeBegin = clamp(newTimeBegin, 0, samples_);
		newTimeEnd = clamp(newTimeEnd, 0, samples_);

		if (newTimeBegin >= timeEnd_ || newTimeEnd <= timeBegin_ ||
			newTimeBegin == newTimeEnd || timeBegin_ == timeEnd_)
		{
			// The time-windows do not share any
			// elements. Clear all objects from
			// the kdtree.

			kdTree_.eraseObjects();
			objectSet_.clear();

			// And insert the new ones.

			kdTree_.insert(
				pointSet_.begin() + newTimeBegin * signals_,
				pointSet_.begin() + newTimeEnd * signals_,
				std::back_inserter(objectSet_));
		}
		else
		{
			// The time-windows overlap. Keep those objects
			// which are common to both, remove and insert
			// appropriately to get the new time-window.

			integer deltaRight = timeEnd_ - newTimeEnd;
			if (deltaRight > 0)
			{
				// Remove objects from the right.

				const integer amount = deltaRight * signals_;
				for (integer i = 0;i < amount;++i)
				{
					kdTree_.erase(objectSet_.back());
					objectSet_.pop_back();
				}
			}
			else if (deltaRight < 0)
			{
				// Add objects to the right.

				kdTree_.insert(
					pointSet_.begin() + timeEnd_ * signals_,
					pointSet_.begin() + newTimeEnd * signals_,
					std::back_inserter(objectSet_));
			}

			integer deltaLeft = timeBegin_ - newTimeBegin;
			if (deltaLeft > 0)
			{
				// Add objects to the left.

				kdTree_.insert(
					boost::make_reverse_iterator(
					pointSet_.begin() + timeBegin_ * signals_),
					boost::make_reverse_iterator(
					pointSet_.begin() + newTimeBegin * signals_),
					std::front_inserter(objectSet_));
			}
			else if (deltaLeft < 0)
			{
				// Remove objects from the left.

				const integer amount = -deltaLeft * signals_;
				for (integer i = 0;i < amount;++i)
				{
					kdTree_.erase(objectSet_.front());
					objectSet_.pop_front();
				}
			}
		}

		timeBegin_ = newTimeBegin;
		timeEnd_ = newTimeEnd;
	}

	const SignalPointSet::KdTree& SignalPointSet::kdTree() const
	{
		return kdTree_;
	}

	SignalPointSet::ConstObjectIterator_Iterator 
		SignalPointSet::begin() const
	{
		return objectSet_.begin();
	}

	SignalPointSet::ConstObjectIterator_Iterator 
		SignalPointSet::end() const
	{
		return objectSet_.end();
	}

	integer SignalPointSet::timeBegin() const
	{
		return timeBegin_;
	}

	integer SignalPointSet::timeEnd() const
	{
		return timeEnd_;
	}

	integer SignalPointSet::samples() const
	{
		return samples_;
	}

	integer SignalPointSet::dimension() const
	{
		return dimension_;
	}

	integer SignalPointSet::dimensionBegin() const
	{
		return dimensionBegin_;
	}

	VectorD SignalPointSet::point(const Object& object) const
	{
		return kdTree_.point(object);
	}

}
