/*
#include "tim/core/mutual_information.h"

#include <pastel/sys/nulliterator.h>

#include <pastel/math/normbijections.h>

using namespace Tim;

namespace
{

	void useCase()
	{
		SignalPtr xSignal;
		SignalPtr ySignal;

		SignalPtr signalSet[2] = {xSignal, ySignal};

		const integer timeWindowRadius = 10;
		const integer kNearest = 1;
		real averageEstimate = 0;

		NullIterator estimateSet;
		Euclidean_NormBijection<real> normBijection;

		// To compute differential entropy for
		// a single signal:

		averageEstimate = mutualInformation(
			xSignal);
		averageEstimate = mutualInformation(
			xSignal);
		averageEstimate = mutualInformation(
			xSignal, kNearest);

		// To compute temporal differential
		// entropy where the neighborhood is given by
		// a time-window:

		temporalMutualInformation(
			xSignal, timeWindowRadius, estimateSet,
			kNearest, normBijection);

		averageEstimate = mutualInformation(
			xSignal, kNearest, normBijection);

		temporalMutualInformation(xSignal, timeWindowRadius, estimateSet, 
			kNearest);
		temporalMutualInformation(xSignal, timeWindowRadius, estimateSet,
			kNearest, normBijection);

		// To compute differential entropy for
		// a set of signals, where each signal is a 
		// different trial of the same experiment:

		averageEstimate = mutualInformation(
			range(signalSet));
		averageEstimate = mutualInformation(
			range(signalSet));
		averageEstimate = mutualInformation(
			range(signalSet), kNearest);
		averageEstimate = mutualInformation(
			range(signalSet) 
			kNearest, normBijection);

		// To compute temporal differential entropy for
		// a set of signals, where each signal is a 
		// different trial of the same experiment:

		temporalMutualInformation(
			range(signalSet), 
			timeWindowRadius, estimateSet); 
		temporalMutualInformation(
			range(signalSet), 
			timeWindowRadius, estimateSet); 
		temporalMutualInformation(
			range(signalSet), 
			timeWindowRadius, estimateSet, 
			kNearest);
		temporalMutualInformation(
			range(signalSet), timeWindowRadius, 
			estimateSet,
			kNearest,
			normBijection);
	}

}
*/