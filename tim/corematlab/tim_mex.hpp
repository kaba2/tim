#ifndef TIM_TIM_MEX_HPP
#define TIM_TIM_MEX_HPP

#include "tim/matlab/tim_mex.h"

#include "tim/core/signal.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/ensure.h>
#include <pastel/sys/sequence/copy_n.h>
#include <pastel/matlab/matlab_argument.h>

namespace Tim
{

	inline Signal asSignal(const mxArray* signal)
	{
		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.

		integer samples = mxGetN(signal);
		integer dimension = mxGetM(signal);
		
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
			if (!isNan(matlabAsScalar<real>(signal, nans * dimension)))
			{
				break;
			}
			++nans;
		}

		if (mxGetClassID(signal) == typeToMatlabClassId<real>())
		{
			// Since the types match, alias the data.
			real* rawData = mxGetPr(signal);
			return Signal(
				samples - nans, 
				dimension, 
				nans, 
				rawData + nans * dimension);
		}

		// The types do not match, so allocate new
		// memory and convert the data to the correct type.
		Signal result(
			samples - nans, 
			dimension, 
			nans);

		matlabGetScalars(signal, result.data().begin(), nans * dimension);

		return result;
	}

	inline std::vector<Signal> getSignals(
		const mxArray* input)
	{
		integer n = 
			mxGetNumberOfElements(input);

		std::vector<Signal> signalSet;
		signalSet.reserve(n);

		for (integer i = 0;i < n;++i)
		{
			const mxArray* signal = mxGetCell(input, i);
			signalSet.push_back(asSignal(signal));
		}

		return signalSet;
	}

	inline Array<Signal> getSignalArray(
		const mxArray* signalSetArray)
	{
		integer signals = mxGetM(signalSetArray);
		integer trials = mxGetN(signalSetArray);
		
		Array<Signal> signalSet(
			Vector2i(trials, signals));

		for (integer y = 0;y < signals;++y)
		{
			for (integer x = 0;x < trials;++x)
			{

				const mxArray* signal = 
					mxGetCell(signalSetArray, signals * x + y);

				signalSet(x, y) = asSignal(signal);
			}
		}

		return signalSet;
	}

}

#endif
