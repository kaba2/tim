#ifndef TIM_SIGNAL_H
#define TIM_SIGNAL_H

#include "tim/core/mytypes.h"

#include <pastel/math/matrix.h>

#include <boost/shared_ptr.hpp>

#include <vector>

namespace Tim
{

	//! A finite sequence in R^n.
	/*!
	A signal is a finite sequence (m-tuple) in R^n. 
	Its dimension is n. Its size is m. The data
	of a signal can be accessed in two ways:
	via a 2-dimensional nxm array view, or
	as a sequence of samples in R^n.
	Both refer to the same memory and thus can
	be used to modify the signal.
	*/

	typedef MatrixD Signal;

	typedef boost::shared_ptr<Signal> SignalPtr;

}

#include "tim/core/signal.hpp"

#endif
