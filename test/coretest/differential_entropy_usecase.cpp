#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/nulliterator.h>

#include <pastel/math/normbijection.h>

using namespace Tim;

namespace
{

	void useCase()
	{
		SignalPtr xSignal = generateGaussian(10000, 3);
		SignalPtr ySignal = generateGaussian(10000, 5);

		SignalPtr signalSet[2] = {xSignal, ySignal};

		const integer timeWindowRadius = 10;
		const real maxRelativeError = 0;
		const integer kNearest = 1;
		real averageEstimate = 0;

		NullIterator estimateSet;
		Euclidean_NormBijection<real> normBijection;

		// To compute differential entropy for
		// a single signal:

		averageEstimate = differentialEntropy(
			xSignal);
		averageEstimate = differentialEntropy(
			xSignal, maxRelativeError);
		averageEstimate = differentialEntropy(
			xSignal, maxRelativeError, kNearest);

		// To compute temporal differential
		// entropy where the neighborhood is given by
		// a time window:

		temporalDifferentialEntropy(
			xSignal, timeWindowRadius, estimateSet,
			maxRelativeError, kNearest, normBijection);

		averageEstimate = differentialEntropy(
			xSignal, maxRelativeError, kNearest, normBijection);

		temporalDifferentialEntropy(xSignal, timeWindowRadius, estimateSet, 
			maxRelativeError, kNearest);
		temporalDifferentialEntropy(xSignal, timeWindowRadius, estimateSet,
			maxRelativeError, kNearest, normBijection);

		// To compute differential entropy for
		// a set of signals, where each signal is a 
		// different trial of the same experiment:

		averageEstimate = differentialEntropy(
			forwardRange(signalSet));
		averageEstimate = differentialEntropy(
			forwardRange(signalSet), maxRelativeError);
		averageEstimate = differentialEntropy(
			forwardRange(signalSet), maxRelativeError, kNearest);
		averageEstimate = differentialEntropy(
			forwardRange(signalSet), maxRelativeError, 
			kNearest, normBijection);

		// To compute temporal differential entropy for
		// a set of signals, where each signal is a 
		// different trial of the same experiment:

		temporalDifferentialEntropy(
			forwardRange(signalSet), 
			timeWindowRadius, estimateSet); 
		temporalDifferentialEntropy(
			forwardRange(signalSet), 
			timeWindowRadius, estimateSet, 
			maxRelativeError); 
		temporalDifferentialEntropy(
			forwardRange(signalSet), 
			timeWindowRadius, estimateSet, 
			maxRelativeError, kNearest);
		temporalDifferentialEntropy(
			forwardRange(signalSet), timeWindowRadius, 
			estimateSet,
			maxRelativeError, kNearest,
			normBijection);
	}

}
