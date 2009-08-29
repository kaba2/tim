#include "tim/core/signalpointset.h"
#include "tim/core/signal_tools.h"

#include <pastel/geometry/search_all_neighbors_pointkdtree.h>
#include <pastel/geometry/slidingmidpoint_splitrule_pointkdtree.h>

#include <pastel/sys/constantiterator.h>

#include <pastel/math/normbijection.h>

namespace Tim
{

	SignalPointSet::SignalPointSet(const SignalPointSet& that)
	{
		ENSURE(false);
	}

	void SignalPointSet::construct(
		integer timeBegin,
		integer timeEnd,
		integer dimensionBegin,
		integer dimensionEnd)
	{
		ENSURE_OP(timeBegin, <=, timeEnd);
		ENSURE_OP(timeBegin, >=, 0);
		ENSURE_OP(timeEnd, <=, samples_);
		ENSURE_OP(dimensionBegin, <=, dimensionEnd);
		ENSURE_OP(dimensionBegin, >=, 0);
		ENSURE_OP(dimensionEnd, <=, signalSet_.front()->dimension());

		dimensionBegin_ = dimensionBegin;
		dimension_ = dimensionEnd - dimensionBegin;

		const integer signals = signalSet_.size();

		// Store points in an interleaved
		// manner.

		pointSet_.reserve(samples_ * signals);

		for (integer t = 0;t < samples_;++t)
		{
			for (integer i = 0;i < signals;++i)
			{
				pointSet_.push_back(
					signalSet_[i]->pointBegin(dimensionBegin_)[t]);
			}
		}

		const integer timeWidth = timeEnd - timeBegin;
		const bool gonnaEraseAll = 
			timeWidth <= samples_ / 2;

		// Insert all the points into the tree.

		if (gonnaEraseAll)
		{
			kdTree_.insert(pointSet_.begin(), pointSet_.end());
		}
		else
		{
			kdTree_.insert(
				pointSet_.begin(), 
				pointSet_.begin() + timeBegin);
			kdTree_.insert(
				pointSet_.begin() + timeBegin, 
				pointSet_.begin() + timeEnd,
				std::back_inserter(objectSet_));
			kdTree_.insert(
				pointSet_.begin() + timeEnd, 
				pointSet_.end());
			
			timeBegin_ = timeBegin;
			timeEnd_ = timeEnd;
		}

		// Compute a fine subdivision for the points.

		kdTree_.refine(SlidingMidpoint2_SplitRule_PointKdTree());

		// Set the initial time window.

		if (gonnaEraseAll)
		{
			kdTree_.eraseObjects();
			setTimeWindow(timeBegin, timeEnd);
		}
	}

	void SignalPointSet::swap(SignalPointSet& that)
	{
		kdTree_.swap(that.kdTree_);
		signalSet_.swap(that.signalSet_);
		pointSet_.swap(that.pointSet_);
		objectSet_.swap(that.objectSet_);
		std::swap(samples_, that.samples_);
		std::swap(timeBegin_, that.timeBegin_);
		std::swap(timeEnd_, that.timeEnd_);
	}

	SignalPointSet& SignalPointSet::operator=(const SignalPointSet& that)
	{
		SignalPointSet copy(that);
		swap(copy);
		return *this;
	}

	void SignalPointSet::setTimeWindow(integer newTimeBegin, integer newTimeEnd)
	{
		ENSURE_OP(newTimeBegin, <=, newTimeEnd);
		ENSURE_OP(newTimeBegin, >=, 0);
		ENSURE_OP(newTimeEnd, <=, samples_);

		if (newTimeBegin == newTimeEnd)
		{
			newTimeBegin = samples_;
			newTimeEnd = samples_;
		}

		const integer signals = signalSet_.size();

		if (newTimeBegin >= timeEnd_ || newTimeEnd <= timeBegin_ ||
			newTimeBegin == newTimeEnd)
		{
			// The time windows do not share any
			// elements. Clear all objects from
			// the kdtree.

			kdTree_.eraseObjects();
			objectSet_.clear();

			// And insert the new ones.

			kdTree_.insert(
				pointSet_.begin() + newTimeBegin * signals,
				pointSet_.begin() + newTimeEnd * signals,
				std::back_inserter(objectSet_));
		}
		else
		{
			// The time windows overlap. Keep those objects
			// which are common to both, remove and insert
			// appropriately to get the new time window.

			integer deltaRight = timeEnd_ - newTimeEnd;
			if (deltaRight > 0)
			{
				// Remove objects from the right.

				const integer amount = deltaRight * signals;
				for (integer i = 0;i < amount;++i)
				{
					kdTree_.erase(objectSet_.back());
					objectSet_.pop_back();
				}
			}
			else if (deltaRight < 0)
			{
				// Adds objects to the right.

				kdTree_.insert(
					pointSet_.begin() + timeEnd_ * signals,
					pointSet_.begin() + newTimeEnd * signals,
					std::back_inserter(objectSet_));
			}

			integer deltaLeft = timeBegin_ - newTimeBegin;
			if (deltaLeft > 0)
			{
				// Add objects to the left.

				kdTree_.insert(
					boost::make_reverse_iterator(
					pointSet_.begin() + timeBegin_ * signals),
					boost::make_reverse_iterator(
					pointSet_.begin() + newTimeBegin * signals),
					std::front_inserter(objectSet_));
			}
			else if (deltaLeft < 0)
			{
				// Remove objects from the left.

				const integer amount = -deltaLeft * signals;
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

	PointD SignalPointSet::point(const Object& object) const
	{
		return kdTree_.point(object);
	}

}
