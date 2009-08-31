/*
#include "tim/core/mutual_information.h"

#include <pastel/sys/nulliterator.h>

#include <pastel/math/normbijection.h>

using namespace Tim;

namespace
{

	void useCase()
	{
		SignalPtr xSignal;
		SignalPtr ySignal;

		SignalPtr signalSet[2] = {xSignal, ySignal};

		const integer timeWindowRadius = 10;
		const real maxRelativeError = 0;
		const integer kNearest = 1;
		real averageEstimate = 0;

		NullIterator estimateSet;
		Euclidean_NormBijection<real> normBijection;

		// To compute differential entropy for
		// a single signal:

		averageEstimate = mutualInformation(
			xSignal);
		averageEstimate = mutualInformation(
			xSignal, maxRelativeError);
		averageEstimate = mutualInformation(
			xSignal, maxRelativeError, kNearest);

		// To compute temporal differential
		// entropy where the neighborhood is given by
		// a time window:

		temporalMutualInformation(
			xSignal, timeWindowRadius, estimateSet,
			maxRelativeError, kNearest, normBijection);

		averageEstimate = mutualInformation(
			xSignal, maxRelativeError, kNearest, normBijection);

		temporalMutualInformation(xSignal, timeWindowRadius, estimateSet, 
			maxRelativeError, kNearest);
		temporalMutualInformation(xSignal, timeWindowRadius, estimateSet,
			maxRelativeError, kNearest, normBijection);

		// To compute differential entropy for
		// a set of signals, where each signal is a 
		// different trial of the same experiment:

		averageEstimate = mutualInformation(
			forwardRange(signalSet));
		averageEstimate = mutualInformation(
			forwardRange(signalSet), maxRelativeError);
		averageEstimate = mutualInformation(
			forwardRange(signalSet), maxRelativeError, kNearest);
		averageEstimate = mutualInformation(
			forwardRange(signalSet), maxRelativeError, 
			kNearest, normBijection);

		// To compute temporal differential entropy for
		// a set of signals, where each signal is a 
		// different trial of the same experiment:

		temporalMutualInformation(
			forwardRange(signalSet), 
			timeWindowRadius, estimateSet); 
		temporalMutualInformation(
			forwardRange(signalSet), 
			timeWindowRadius, estimateSet, 
			maxRelativeError); 
		temporalMutualInformation(
			forwardRange(signalSet), 
			timeWindowRadius, estimateSet, 
			maxRelativeError, kNearest);
		temporalMutualInformation(
			forwardRange(signalSet), timeWindowRadius, 
			estimateSet,
			maxRelativeError, kNearest,
			normBijection);
	}

}
*/