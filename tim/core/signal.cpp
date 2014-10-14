#include "tim/core/signal.h"

#include <pastel/sys/subview.h>
#include <pastel/sys/view_tools.h>

namespace Tim
{

	Signal::Signal()
		: data_(0, 0)
		, t_(0)
	{
	}

	Signal::Signal(const Signal& that)
	: data_(that.data_)
	, t_(that.t_)
	{
	}

	Signal::Signal(Signal&& that)
	: Signal()
	{
		swap(that);
	}

	Signal::Signal(
		integer samples, integer dimension, integer t)
		: data_(samples, dimension)
		, t_(t)
	{
	}

	Signal::Signal(
		integer samples, integer dimension, 
		integer t, real* dataToAlias)
		: data_(samples, dimension, withAliasing(dataToAlias))
		, t_(t)
	{
	}

	void Signal::swap(Signal& that)
	{
		data_.swap(that.data_);
		std::swap(t_, that.t_);
	}

	Signal& Signal::operator=(Signal that)
	{
		swap(that);
		return *this;
	}

	integer Signal::dimension() const
	{
		return data_.width();
	}

	integer Signal::samples() const
	{
		return data_.height();
	}

	Matrix<real>& Signal::data()
	{
		return data_;
	}

	const Matrix<real>& Signal::data() const
	{
		return data_;
	}

	void Signal::setT(integer t)
	{
		t_ = t;
	}

	integer Signal::t() const
	{
		return t_;
	}

	Signal::Point_Iterator Signal::pointBegin(integer dimensionBegin)
	{
		return sparseIterator(countingIterator(
			data_.rawBegin() + dimensionBegin), dimension());
	}

	Signal::Point_ConstIterator Signal::pointBegin(integer dimensionBegin) const
	{
		return constSparseIterator(countingIterator(
			data_.rawBegin() + dimensionBegin), dimension());
	}

	Signal::Point_Iterator Signal::pointEnd(integer dimensionBegin)
	{
		return sparseIterator(countingIterator(
			data_.rawBegin() + dimensionBegin + samples() * dimension()), 
			dimension());
	}

	Signal::Point_ConstIterator Signal::pointEnd(integer dimensionBegin) const
	{
		return constSparseIterator(countingIterator(
			data_.rawBegin() + dimensionBegin + samples() * dimension()), 
			dimension());
	}

}
