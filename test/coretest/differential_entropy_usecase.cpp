#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/iterator/null_iterator.h>

#include <pastel/math/normbijection.h>

using namespace Tim;

namespace
{

	void useCase()
	{
		SignalData xSignal = generateGaussian(3, 10000);
		SignalData ySignal = generateGaussian(5, 10000);

		Signal signalSet[2] = {xSignal, ySignal};

		const integer timeWindowRadius = 10;
		const integer kNearest = 1;
		dreal averageEstimate = 0;

		Euclidean_Norm<dreal> norm;

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
			kNearest, norm);

		averageEstimate = differentialEntropyKl(
			constantRange(xSignal), kNearest, norm);

		temporalDifferentialEntropyKl(
			constantRange(xSignal), timeWindowRadius,
			kNearest);
		temporalDifferentialEntropyKl(
			constantRange(xSignal), timeWindowRadius,
			kNearest, norm);

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
			kNearest, norm);

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
			norm);
	}

}
