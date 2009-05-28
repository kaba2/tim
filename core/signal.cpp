#include "tim/core/signal.h"

#include <pastel/sys/subview.h>
#include <pastel/sys/view_tools.h>

namespace Tim
{

	Signal::Signal(integer dimension, integer size)
		: data_(size, dimension, 0)
		, pointSet_()
	{
		ENSURE1(dimension >= 0, dimension);
		ENSURE1(size >= 0, size);

		updatePointSet();
	}

	Signal::Signal(
		const SignalPtr& signalToAlias,
		integer dimensionBegin,
		integer dimensionEnd)
		: data_(signalToAlias->size(), std::abs(dimensionEnd - dimensionBegin), 
		withAliasing(&signalToAlias->view()(0, dimensionBegin)))
		, pointSet_()
	{
		// Note the mabs() above is just because we want 
		// to catch precondition violations here.

		ENSURE2(dimensionBegin >= 0 && dimensionBegin < signalToAlias->dimension(), 
			dimensionBegin, signalToAlias->dimension());
		ENSURE2(dimensionEnd >= 0 && dimensionEnd <= signalToAlias->dimension(), 
			dimensionEnd, signalToAlias->dimension());
		ENSURE2(dimensionBegin <= dimensionEnd, dimensionBegin, dimensionEnd);

		updatePointSet();
	}

	void Signal::swap(Signal& that)
	{
		data_.swap(that.data_);
		pointSet_.swap(that.pointSet_);
	}

	integer Signal::size() const
	{
		return data_.height();
	}

	integer Signal::dimension() const
	{
		return data_.width();
	}

	SignalView Signal::view()
	{
		return arrayView(data_);
	}

	ConstSignalView Signal::constView() const
	{
		return constArrayView(data_);
	}

	const std::vector<DynamicPoint>& Signal::pointSet() const
	{
		return pointSet_;
	}

	DynamicPoint& Signal::operator[](integer index)
	{
		PENSURE2(index >= 0 && index < pointSet_.size(), index, pointSet_.size());

		return pointSet_[index];
	}

	const DynamicPoint& Signal::operator[](integer index) const
	{
		PENSURE2(index >= 0 && index < pointSet_.size(), index, pointSet_.size());

		return pointSet_[index];
	}

	// Private

	void Signal::updatePointSet()
	{
		const integer size = data_.width();
		const integer dimension = data_.height();

		pointSet_.resize(size, nullPoint<Dynamic, real>());

		for (integer i = 0;i < size;++i)
		{
			DynamicPoint temp(
				ofDimension(dimension),
				withAliasing(&data_(i, 0)));

			pointSet_[i].swap(temp);

			ASSERT(&pointSet_[i][0] == &data_(0, i));
		}
	}

}
