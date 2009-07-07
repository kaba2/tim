#ifndef TIM_PARTIAL_MUTUAL_INFORMATION_H
#define TIM_PARTIAL_MUTUAL_INFORMATION_H

#include "tim/core/signal.h"

namespace Tim
{

	//! Computes partial mutual information I(A, B | C).
	/*!
	Preconditions:
	kNearest > 0
	aSignal->samples() == bSignal->samples() 
	aSignal->samples() == cSignal->samples()

	See:
	"Partial Mutual Information for Coupling Analysis 
	of Multivariate Time Series", 
	Stefan Frenzel and Bernd Pompe,
	Physical Review Letters, 2007.
	*/

	TIMCORE real partialMutualInformation(
		const SignalPtr& aSignal,
		const SignalPtr& bSignal,
		const SignalPtr& cSignal,
		integer kNearest);

}

#endif
