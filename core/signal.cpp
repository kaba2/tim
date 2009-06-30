#include "tim/core/signal.h"

#include <pastel/sys/subview.h>
#include <pastel/sys/view_tools.h>

namespace Tim
{

	/*

	Signal::Signal(integer dimension, integer size, real defaultData)
		: data_(size, dimension)
	{
		ENSURE1(dimension >= 0, dimension);
		ENSURE1(size >= 0, size);

		data_.set(defaultData);
	}

	void Signal::swap(Signal& that)
	{
		data_.swap(that.data_);
	}

	integer Signal::size() const
	{
		return data_.height();
	}

	integer Signal::width() const
	{
		return data_.width();
	}

	SignalView Signal::view()
	{
		return data_.view();
	}

	ConstSignalView Signal::constView() const
	{
		return data_.constView();
	}

	MatrixD::Row Signal::operator[](integer index)
	{
		return data_[index];
	}

	MatrixD::ConstRow Signal::operator[](integer index) const
	{
		return data_[index];
	}
	*/

}
