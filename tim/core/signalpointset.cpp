#include "tim/core/signalpointset.h"
#include "tim/core/signal_tools.h"

#include <pastel/geometry/slidingmidpoint_splitrule.h>
#include <pastel/geometry/difference_alignedbox_alignedbox.h>
#include <pastel/geometry/intersect_alignedbox_alignedbox.h>
#include <pastel/geometry/overlaps_alignedbox_alignedbox.h>
#include <pastel/geometry/contains_alignedbox_alignedbox.h>

#include <pastel/sys/constant_iterator.h>

#include <pastel/math/normbijections.h>

#include <boost/iterator/reverse_iterator.hpp>
#include <boost/bind.hpp>

namespace Tim
{

	SignalPointSet::SignalPointSet(const SignalPointSet& that)
	{
		ENSURE(false);
	}

	void SignalPointSet::swap(SignalPointSet& that)
	{
		kdTree_.swap(that.kdTree_);
		pointSet_.swap(that.pointSet_);
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

		const AlignedBox<integer, 1> sampleWindow(
			timeBegin_, timeBegin_ + samples_);
		const AlignedBox<integer, 1> window(
			windowBegin_, windowEnd_);
		AlignedBox<integer, 1> newWindow(
			newWindowBegin, newWindowEnd);

		// Cut the new window to the defined time range.
		if (!intersect(newWindow, sampleWindow, newWindow))
		{
			// The new window does not overlap with the
			// sample window. Hide all points.
			kdTree_.hide();
		}
		else
		{
			if (!overlaps(window, newWindow))
			{
				// The new window does not contain any of the
				// existing points.
				kdTree_.hide();
			}
			else
			{			
				// Hide those points which are not in the new window.
				difference(window, newWindow, 
					boost::bind(&SignalPointSet::hide, this, _1));
			}

			if (contains(newWindow, sampleWindow))
			{
				// The new window contains all points.
				kdTree_.show();
			}
			else
			{
				// Show all those points not yet in the window.
				difference(newWindow, window, 
					boost::bind(&SignalPointSet::show, this, _1));
			}
		}

		windowBegin_ = newWindow.min().x();
		windowEnd_ = newWindow.max().x();
	}

	const SignalPointSet::KdTree& SignalPointSet::kdTree() const
	{
		return kdTree_;
	}

	SignalPointSet::Point_ConstIterator_Iterator 
		SignalPointSet::begin() const
	{
		return pointSet_.begin() + (windowBegin_ - timeBegin_) * signals_;
	}

	SignalPointSet::Point_ConstIterator_Iterator 
		SignalPointSet::end() const
	{
		return pointSet_.begin() + (windowEnd_ - timeBegin_) * signals_;
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
		return pointAsVector(point, kdTree_.locator());
	}

	// Private

	void SignalPointSet::hide(
		const AlignedBox<integer, 1>& range)
	{
		const integer iBegin = (range.min().x() - timeBegin_) * signals_;
		const integer iEnd = (range.max().x() - timeBegin_) * signals_;
		for (integer i = iBegin;i < iEnd;++i)
		{
			kdTree_.hide(pointSet_[i]);
		}
	}

	void SignalPointSet::show(
		const AlignedBox<integer, 1>& range)
	{
		const integer iBegin = (range.min().x() - timeBegin_) * signals_;
		const integer iEnd = (range.max().x() - timeBegin_) * signals_;
		for (integer i = iBegin;i < iEnd;++i)
		{
			kdTree_.show(pointSet_[i]);
		}
	}

}
