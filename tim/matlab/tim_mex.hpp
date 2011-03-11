#ifndef TIM_TIM_MEX_HPP
#define TIM_TIM_MEX_HPP

#include "tim/matlab/tim_mex.h"

#include "tim/core/signal.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/ensure.h>
#include <pastel/sys/copy_n.h>
#include <pastel/matlab/matlab_argument.h>

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
		
		/*!
		The output signal aliases the input signal such that
		the NaNs from the beginning of the signal are skipped.

		TIM Matlab allows two forms for signals, a _lagged form_, 
		and a _NaN-padded form_. In NaN-padded form the start of the 
		signal contains NaNs, whose role is only to keep the first
		sample associated to the time instant zero. When converting 
		from NaN-padded form to lagged form, the NaN part is removed 
		from the signal and the signal is interpreted to begin at
		the first non-NaN sample. To maintain time-correspondence,
		the first sample in the lagged form is taken to begin at a 
		time instant given by the number of NaN samples in the 
		NaN-padded form.
		*/

		integer nans = 0;
		while(nans < samples)
		{
			if (!isNan(asScalar<real>(signal, nans * dimension)))
			{
				break;
			}
			++nans;
		}

		SignalPtr result;
		if (mxGetClassID(signal) == typeToMatlabClassId<real>())
		{
			// Since the types match, alias the data.
			real* rawData = mxGetPr(signal);
			result = SignalPtr(
				new Signal(samples - nans, dimension, 
				nans, rawData + nans * dimension));
		}
		else
		{
			// The types do not match, so allocate new
			// memory and convert the data to the correct type.
			result = SignalPtr(
				new Signal(samples - nans, dimension, nans));
			getScalars(signal, result->data().begin(), nans * dimension);
		}

		return result;
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
