// Description: Mutual information estimation via Kraskow's algorithm

#ifndef TIM_MUTUAL_INFORMATION_KRASKOW_H
#define TIM_MUTUAL_INFORMATION_KRASKOW_H

#include "tim/core/signal.h"

#include <pastel/math/matrix.h>
#include <pastel/math/cholesky_decomposition.h>

#include <pastel/sys/smallset.h>

namespace Tim
{

	//! Computes mutual information between a set of signals.
	/*!
	Preconditions:
	kNearest > 0
	maxRelativeError >= 0
	Signal_ForwardIterator dereferences to SignalPtr.
	
	signalBegin, signalEnd:
	A set of signals. If the number of samples varies between
	signals, then the minimum number of samples among the signals
	is used to compute the mi.

	kNearest:
	The number of nearest neighbors to use in the estimation.

	maxRelativeError:
	The maximum relative error in the distance of the
	nearest neighbors. Allowing approximate results in
	nearest neighbor searching can accelerate performance 
	by one to two orders of magnitude in higher dimensions.
	0 = exact, 1 = 100% relative error, 2 = 200% relative error,
	etc.

	By 'mutual information' we actually mean 'total correlation' 
	which is one generalization of mutual information for higher 
	number of signals than two.
	If H(X) denotes the differential entropy of a real 
	random variable X, and we are given n random variables
	(X_1, ..., X_n)	then the total correlation between 
	them is defined by:
	I(X_1, ..., X_n) = (sum_i H(X_i)) - H(X_1, ..., X_n)
	While differential entropy does not make
	sense w.r.t the amount of information, total correlation does.
	See 'tim/core/differential_entropy.h' for more information
	on differential entropy.
	*/

	template <
		typename Signal_A_ForwardIterator,
		typename Signal_B_ForwardIterator,
		typename Lag_ForwardIterator>
	real mutualInformation(
		const Signal_A_ForwardIterator& aTrialBegin,
		const Signal_B_ForwardIterator& bTrialBegin,
		integer trials,
		const Lag_ForwardIterator& lagBegin,
		integer lags,
		integer sigma = -1,
		integer kNearest = 1)

}

#include "tim/core/mutual_information_kraskow.hpp"

#endif
