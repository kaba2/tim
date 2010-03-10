#ifndef TIM_TIM_MEX_HPP
#define TIM_TIM_MEX_HPP

#include "tim/matlab/tim_mex.h"

#include "tim/core/signal.h"

#include <pastel/sys/ensure.h>
#include <pastel/sys/pastelomp.h>
#include <pastel/sys/stdext_copy_n.h>

namespace Tim
{

	template <typename SignalPtr_OutputIterator>
	void getSignals(const mxArray* input,
					SignalPtr_OutputIterator output)
	{
		const integer trials = mxGetNumberOfElements(input);

		for (integer i = 0;i < trials;++i)
		{
			const mxArray* signalArray = mxGetCell(input, i);

			// It is intentional to assign the width
			// and height the wrong way. The reason
			// is that Matlab uses column-major storage
			// while we use row-major storage.
			const integer samples = mxGetN(signalArray);
			const integer dimension = mxGetM(signalArray);

			real* rawData = mxGetPr(signalArray);

			*output = SignalPtr(
				new Signal(samples, dimension, rawData));
			++output;
		}
	}

	template <typename Integer_OutputIterator>
	void getIntegers(const mxArray* input,
					 Integer_OutputIterator output)
	{
		getReals(input, output);
	}

	template <typename Integer_OutputIterator>
	void getReals(const mxArray* input,
				  Integer_OutputIterator output)
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
		Array<SignalPtr, 2>& signalSet)
	{
		const integer signals = mxGetM(signalSetArray);
		const integer trials = mxGetN(signalSetArray);
		
		signalSet.setExtent(trials, signals);

		for (integer y = 0;y < signals;++y)
		{
			for (integer x = 0;x < trials;++x)
			{
				const mxArray* signalArray = mxGetCell(signalSetArray, signals * x + y);

				const integer samples = mxGetN(signalArray);
				const integer dimension = mxGetM(signalArray);

				real* rawData = mxGetPr(signalArray);

				signalSet(x, y) = SignalPtr(
					new Signal(samples, dimension, rawData));
			}
		}
	}

}

#endif
