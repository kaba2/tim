#ifndef TIM_SIGNAL_H
#define TIM_SIGNAL_H

#include "tim/core/mytypes.h"

#include <pastel/math/matrix.h>

#include <pastel/sys/countedptr.h>

#include <vector>

namespace Tim
{

	//! A time series.

	class TIMCORE Signal
		: public ReferenceCounted
	{
	public:
		// Using default copy constructor.
		// Using default assignment.
		// Using default destructor.

		Signal();
		Signal(integer samples, integer dimension);
		Signal(integer samples, integer dimension,
			real* dataToAlias);

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

	typedef CountedPtr<Signal> SignalPtr;

}

#include "tim/core/signal.hpp"

#endif
