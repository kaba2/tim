#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/null_iterator.h>

#include <pastel/math/normbijections.h>

using namespace Tim;

namespace
{

	void useCase()
	{
		SignalPtr xSignal = generateGaussian(10000, 3);
		SignalPtr ySignal = generateGaussian(10000, 5);

		SignalPtr signalSet[2] = {xSignal, ySignal};

		const integer timeWindowRadius = 10;
		const integer kNearest = 1;
		real averageEstimate = 0;

		Euclidean_NormBijection<real> normBijection;

		// To compute differential entropy for
		// a single signal:

		averageEstimate = differentialEntropyKl(
			constantRange(xSignal));
		averageEstimate = differentialEntropyKl(
			constantRange(xSignal));
		averageEstimate = differentialEntropyKl(
			constantRange(xSignal), kNearest);

		// To compute temporal differential
		// entropy where the neighborhood is given by
		// a time-window:

		temporalDifferentialEntropyKl(
			constantRange(xSignal), timeWindowRadius,
			kNearest, normBijection);

		averageEstimate = differentialEntropyKl(
			constantRange(xSignal), kNearest, normBijection);

		temporalDifferentialEntropyKl(
			constantRange(xSignal), timeWindowRadius,
			kNearest);
		temporalDifferentialEntropyKl(
			constantRange(xSignal), timeWindowRadius,
			kNearest, normBijection);

		// To compute differential entropy for
		// a set of signals, where each signal is a 
		// different trial of the same experiment:

		averageEstimate = differentialEntropyKl(
			range(signalSet));
		averageEstimate = differentialEntropyKl(
			range(signalSet));
		averageEstimate = differentialEntropyKl(
			range(signalSet), kNearest);
		averageEstimate = differentialEntropyKl(
			range(signalSet), 
			kNearest, normBijection);

		// To compute temporal differential entropy for
		// a set of signals, where each signal is a 
		// different trial of the same experiment:

		temporalDifferentialEntropyKl(
			range(signalSet), 
			timeWindowRadius); 
		temporalDifferentialEntropyKl(
			range(signalSet), 
			timeWindowRadius); 
		temporalDifferentialEntropyKl(
			range(signalSet), 
			timeWindowRadius, 
			kNearest);
		temporalDifferentialEntropyKl(
			range(signalSet), timeWindowRadius, 
			kNearest,
			normBijection);
	}

}
