#include "tim/core/signal.h"

#include <pastel/sys/subview.h>
#include <pastel/sys/view_tools.h>

namespace Tim
{

	Signal::Signal()
		: name_()
		, data_()
	{
	}

	Signal::Signal(integer samples, integer dimension)
		: name_()
		, data_(samples, dimension)
	{
	}

	Signal::Signal(integer samples, integer dimension, 
		real* dataToAlias)
		: name_()
		, data_(samples, dimension, withAliasing(dataToAlias))
	{
	}

	void Signal::setName(const std::string& name)
	{
		name_ = name;
	}

	const std::string& Signal::name() const
	{
		return name_;
	}

	integer Signal::dimension() const
	{
		return data_.width();
	}

	integer Signal::samples() const
	{
		return data_.height();
	}

	MatrixD& Signal::data()
	{
		return data_;
	}

	const MatrixD& Signal::data() const
	{
		return data_;
	}

	Signal::PointIterator Signal::pointBegin(integer dimensionBegin)
	{
		return sparseIterator(countingIterator(
			data_.rawBegin() + dimensionBegin), dimension());
	}

	Signal::ConstPointIterator Signal::pointBegin(integer dimensionBegin) const
	{
		return constSparseIterator(countingIterator(
			data_.rawBegin() + dimensionBegin), dimension());
	}

	Signal::PointIterator Signal::pointEnd(integer dimensionBegin)
	{
		return sparseIterator(countingIterator(
			data_.rawBegin() + dimensionBegin + data_.size()), 
			dimension());
	}

	Signal::ConstPointIterator Signal::pointEnd(integer dimensionBegin) const
	{
		return constSparseIterator(countingIterator(
			data_.rawBegin() + dimensionBegin + data_.size()), 
			dimension());
	}

}
