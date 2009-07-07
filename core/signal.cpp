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

	Signal::Signal(integer dimension, integer samples,
		real* dataToAlias)
		: name_()
		, data_(dimension, samples, withAliasing(dataToAlias))
	{
	}

	Signal::Signal(integer dimension, integer samples)
		: name_()
		, data_(dimension, samples)
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
		return data_.height();
	}

	integer Signal::samples() const
	{
		return data_.width();
	}

	MatrixD& Signal::data()
	{
		return data_;
	}

	const MatrixD& Signal::data() const
	{
		return data_;
	}

}
