#include "tim/core/signalpointset.h"
#include "tim/core/signal_tools.h"

#include <pastel/geometry/search_all_neighbors_pointkdtree.h>
#include <pastel/geometry/slidingmidpoint_splitrule_pointkdtree.h>

#include <pastel/sys/constantiterator.h>

#include <pastel/math/normbijections.h>

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
			kdTree_.insert(
				range(pointSet_.begin(), pointSet_.end()), 
				std::back_inserter(activeSet_));
			
			windowBegin_ = timeBegin_;
			windowEnd_ = timeBegin_ + pointSet_.size();
		}
		else
		{
			kdTree_.insert(
				range(pointSet_.begin(), pointSet_.end()));
		}

		// Compute a fine subdivision for the points.

		kdTree_.refine(SplitRule());

		if (!startFull)
		{
			kdTree_.erase();
		}
	}

	void SignalPointSet::swap(SignalPointSet& that)
	{
		kdTree_.swap(that.kdTree_);
		pointSet_.swap(that.pointSet_);
		activeSet_.swap(that.activeSet_);
		std::swap(signals_, that.signals_);
		std::swap(samples_, that.samples_);
		std::swap(windowBegin_, that.windowBegin_);
		std::swap(windowEnd_, that.windowEnd_);
		std::swap(dimensionBegin_, that.dimensionBegin_);
		std::swap(dimension_, that.dimension_);
		std::swap(timeBegin_, that.timeBegin_);
	}

	SignalPointSet& SignalPointSet::operator=(
		const SignalPointSet& that)
	{
		SignalPointSet copy(that);
		swap(copy);
		return *this;
	}

	void SignalPointSet::setTimeWindow(
		integer newWindowBegin, integer newWindowEnd)
	{
		ENSURE_OP(newWindowBegin, <=, newWindowEnd);

		newWindowBegin = clamp(newWindowBegin, timeBegin_, timeBegin_ + samples_);
		newWindowEnd = clamp(newWindowEnd, timeBegin_, timeBegin_ + samples_);

		if (newWindowBegin >= windowEnd_ || newWindowEnd <= windowBegin_ ||
			newWindowBegin == newWindowEnd || windowBegin_ == windowEnd_)
		{
			// The time-windows do not share any
			// elements. Clear all points from
			// the kdtree.

			kdTree_.erase();
			activeSet_.clear();

			// And insert the new ones.

			kdTree_.insert(
				range(
				pointSet_.begin() + (newWindowBegin - timeBegin_) * signals_,
				pointSet_.begin() + (newWindowEnd - timeBegin_) * signals_),
				std::back_inserter(activeSet_));
		}
		else
		{
			// The time-windows overlap. Keep those points
			// which are common to both, remove and insert
			// appropriately to get the new time-window.

			integer deltaRight = windowEnd_ - newWindowEnd;
			if (deltaRight > 0)
			{
				// Remove points from the right.

				const integer amount = deltaRight * signals_;
				for (integer i = 0;i < amount;++i)
				{
					kdTree_.erase(activeSet_.back());
					activeSet_.pop_back();
				}
			}
			else if (deltaRight < 0)
			{
				// Add points to the right.

				kdTree_.insert(
					range(
					pointSet_.begin() + (windowEnd_ - timeBegin_) * signals_,
					pointSet_.begin() + (newWindowEnd - timeBegin_) * signals_),
					std::back_inserter(activeSet_));
			}

			integer deltaLeft = windowBegin_ - newWindowBegin;
			if (deltaLeft > 0)
			{
				// Add points to the left.

				kdTree_.insert(
					range(
					boost::make_reverse_iterator(
					pointSet_.begin() + (windowBegin_ - timeBegin_) * signals_),
					boost::make_reverse_iterator(
					pointSet_.begin() + (newWindowBegin - timeBegin_) * signals_)),
					std::front_inserter(activeSet_));
			}
			else if (deltaLeft < 0)
			{
				// Remove points from the left.

				const integer amount = -deltaLeft * signals_;
				for (integer i = 0;i < amount;++i)
				{
					kdTree_.erase(activeSet_.front());
					activeSet_.pop_front();
				}
			}
		}

		windowBegin_ = newWindowBegin;
		windowEnd_ = newWindowEnd;
	}

	const SignalPointSet::KdTree& SignalPointSet::kdTree() const
	{
		return kdTree_;
	}

	SignalPointSet::Point_ConstIterator_Iterator 
		SignalPointSet::begin() const
	{
		return activeSet_.begin();
	}

	SignalPointSet::Point_ConstIterator_Iterator 
		SignalPointSet::end() const
	{
		return activeSet_.end();
	}

	integer SignalPointSet::windowBegin() const
	{
		return windowBegin_;
	}

	integer SignalPointSet::windowEnd() const
	{
		return windowEnd_;
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

	VectorD SignalPointSet::point(const Point& point) const
	{
		return kdTree_.pointPolicy()(point);
	}

}
