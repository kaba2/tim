// Description: Algorithms for Signal's

#ifndef TIM_SIGNAL_TOOLS_H
#define TIM_SIGNAL_TOOLS_H

#include "tim/core/signal.h"
#include "tim/core/signal_merge.h"
#include "tim/core/signal_split.h"
#include "tim/core/signal_properties.h"

#include <pastel/math/matrix/matrix.h>

#include <pastel/sys/range.h>
#include <pastel/sys/array/array.h>

#include <iostream>

namespace Tim
{

	//! Prints the signal into an output stream.
	TIM std::ostream& operator<<(
		std::ostream& stream, const Signal& signal);

	//! Computes the covariance of the signal samples.
	TIM void computeCovariance(
		const Signal& signal,
		Matrix<real>& result);

	//! Transforms the given signal to identity covariance.
	TIM void normalizeCovariance(
		Signal& signal,
		const Matrix<real>& covariance);

}

#include "tim/core/signal_generate.h"

#endif
