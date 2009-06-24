#include "tim/core/signal.h"

#include <pastel/sys/subview.h>
#include <pastel/sys/view_tools.h>

namespace Tim
{

	Signal::Signal(integer dimension, integer size, real defaultData)
		: data_(dimension, size, defaultData)
		, pointSet_()
		, dimensionBegin_(0)
		, dimension_(dimension)
	{
		ENSURE1(dimension >= 0, dimension);
		ENSURE1(size >= 0, size);

		updatePointSet();
	}

	Signal::Signal(
		const SignalPtr& signalToAlias,
		integer dimension,
		integer dimensionBegin)
		: data_(signalToAlias->dimension(), signalToAlias->size(),
		withAliasing(&signalToAlias->view()(0, 0)))
		, pointSet_()
		, dimensionBegin_(dimensionBegin)
		, dimension_(dimension)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(dimensionBegin >= 0, dimensionBegin);
		ENSURE2(dimensionBegin + dimension <= signalToAlias->dimension(),
			dimensionBegin + dimension, signalToAlias->dimension());

		updatePointSet();
	}

	void Signal::swap(Signal& that)
	{
		data_.swap(that.data_);
		pointSet_.swap(that.pointSet_);
		std::swap(dimensionBegin_, that.dimensionBegin_);
		std::swap(dimension_, that.dimension_);
	}

	integer Signal::size() const
	{
		return data_.height();
	}

	integer Signal::dimension() const
	{
		return dimension_;
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
		const integer points = size();

		pointSet_.resize(points, nullPoint<Dynamic, real>());

		for (integer i = 0;i < points;++i)
		{
			DynamicPoint temp(
				ofDimension(dimension_),
				withAliasing(&data_(dimensionBegin_, i)));

			pointSet_[i].swap(temp);

			ASSERT(&pointSet_[i][0] == &data_(dimensionBegin_, i));
		}
	}

}
