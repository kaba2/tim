#ifndef TIM_TIM_MEX_HPP
#define TIM_TIM_MEX_HPP

#include "tim/matlab/tim_mex.h"

#include "tim/core/signal.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/ensure.h>
#include <pastel/sys/copy_n.h>

namespace Tim
{

	inline SignalPtr asSignal(const mxArray* signal)
	{
		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.
		const integer samples = mxGetN(signal);
		const integer dimension = mxGetM(signal);

		real* rawData = mxGetPr(signal);

		const SignalPtr nanForm = SignalPtr(
			new Signal(samples, dimension, 0, rawData));

		return nanToLagged(nanForm);
	}

	template <typename SignalPtr_Iterator>
	integer getSignals(const mxArray* input,
					SignalPtr_Iterator output)
	{
		const integer trials = 
			mxGetNumberOfElements(input);

		for (integer i = 0;i < trials;++i)
		{
			const mxArray* signal = mxGetCell(input, i);
			*output = asSignal(signal);
			++output;
		}

		return trials;
	}

	inline void getSignalArray(
		const mxArray* signalSetArray, 
		Array<SignalPtr>& signalSet)
	{
		const integer signals = mxGetM(signalSetArray);
		const integer trials = mxGetN(signalSetArray);
		
		signalSet.setExtent(trials, signals);

		for (integer y = 0;y < signals;++y)
		{
			for (integer x = 0;x < trials;++x)
			{
				const mxArray* signal = 
					mxGetCell(signalSetArray, signals * x + y);

				signalSet(x, y) = asSignal(signal);
			}
		}
	}

}

#endif
