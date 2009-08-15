#include "tim/core/signalpointset.h"

#include <pastel/geometry/search_all_neighbors_pointkdtree.h>

#include <pastel/sys/constantiterator.h>

#include <pastel/math/normbijection.h>

namespace Tim
{

	SignalPointSet::SignalPointSet(const SignalPointSet& that)
	{
		ENSURE(false);
	}

	void SignalPointSet::swap(SignalPointSet& that)
	{
		kdTree_.swap(that.kdTree_);
		signalSet_.swap(that.signalSet_);
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

			for (integer i = 0;i < signals;++i)
			{
				kdTree_.insert(
					signalSet_[i]->pointBegin() + newTimeBegin,
					signalSet_[i]->pointBegin() + newTimeEnd,
					std::back_inserter(objectSet_));
			}
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

				for (integer i = 0;i < signals;++i)
				{
					kdTree_.insert(
						signalSet_[i]->pointBegin() + timeEnd_,
						signalSet_[i]->pointBegin() + newTimeEnd,
						std::back_inserter(objectSet_));
				}
			}

			integer deltaLeft = timeBegin_ - newTimeBegin;
			if (deltaLeft > 0)
			{
				// Add objects to the left.

				for (integer i = 0;i < signals;++i)
				{
					kdTree_.insert(
						boost::make_reverse_iterator(signalSet_[i]->pointBegin() + timeBegin_),
						boost::make_reverse_iterator(signalSet_[i]->pointBegin() + newTimeBegin),
						std::front_inserter(objectSet_));
				}
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

}
