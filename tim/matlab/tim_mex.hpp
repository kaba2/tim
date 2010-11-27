#ifndef TIM_TIM_MEX_HPP
#define TIM_TIM_MEX_HPP

#include "tim/matlab/tim_mex.h"

#include "tim/core/signal.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/ensure.h>
#include <pastel/sys/pastelomp.h>
#include <pastel/sys/stdext_copy_n.h>

namespace Tim
{

	inline SignalPtr getSignal(const mxArray* signal)
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

	template <typename SignalPtr_OutputIterator>
	void getSignals(const mxArray* input,
					SignalPtr_OutputIterator output)
	{
		const integer trials = mxGetNumberOfElements(input);

		for (integer i = 0;i < trials;++i)
		{
			const mxArray* signal = mxGetCell(input, i);
			*output = getSignal(signal);
			++output;
		}
	}

	template <typename Integer_OutputIterator>
	void getIntegers(const mxArray* input,
					 Integer_OutputIterator output)
	{
		getReals(input, output);
	}

	template <typename Real_OutputIterator>
	void getReals(const mxArray* input,
				  Real_OutputIterator output)
	{
		const integer n = mxGetNumberOfElements(input);
		StdExt::copy_n(mxGetPr(input), n, output);
	}

	inline integer getInteger(const mxArray* input)
	{
		return *mxGetPr(input);
	}

	inline real getReal(const mxArray* input)
	{
		return *mxGetPr(input);
	}

	inline std::string getString(const mxArray* input)
	{
		char* text = mxArrayToString(input);
		std::string result(text);
		mxFree(text);
		return result;
	}

	inline void setNumberOfThreads(integer threads)
	{
		ENSURE_OP(threads, >, 0);
#if PASTEL_ENABLE_OMP != 0
		omp_set_num_threads(threads);
#endif
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

				signalSet(x, y) = getSignal(signal);
			}
		}
	}

}

#endif
