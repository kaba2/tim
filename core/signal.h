#ifndef TIM_SIGNAL_H
#define TIM_SIGNAL_H

#include "tim/core/mytypes.h"

#include <pastel/math/matrix.h>

#include <boost/shared_ptr.hpp>

#include <vector>

namespace Tim
{

	//! A time series.

	class TIMCORE Signal
	{
	public:
		// Using default copy constructor.
		// Using default assignment.
		// Using default destructor.

		Signal();
		Signal(integer dimension, integer samples,
			real* dataToAlias);
		Signal(integer dimension, integer samples);

		void setName(const std::string& name);
		const std::string& name() const;

		integer dimension() const;
		integer samples() const;

		MatrixD& data();
		const MatrixD& data() const;

	private:
		std::string name_;
		MatrixD data_;
	};

	typedef boost::shared_ptr<Signal> SignalPtr;

}

#include "tim/core/signal.hpp"

#endif
