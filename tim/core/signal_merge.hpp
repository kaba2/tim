#ifndef TIM_SIGNAL_MERGE_HPP
#define TIM_SIGNAL_MERGE_HPP

#include "tim/core/signal_merge.h"
#include "tim/core/signal_properties.h"

namespace Tim
{

	template <
		typename SignalPtr_Iterator,
		typename Integer_Iterator>
	SignalPtr merge(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const ForwardRange<Integer_Iterator>& lagSet,
		integer* tMin)
	{
		ENSURE_OP(signalSet.size(), ==, lagSet.size());

		if (tMin)
		{
			*tMin = 0;
		}

		if (signalSet.empty() ||
			lagSet.empty())
		{
			return SignalPtr();
		}

		// Compute joint dimension.

		integer jointDimension = 0;
		SignalPtr_Iterator signalIter = signalSet.begin();
		const SignalPtr_Iterator signalIterEnd = signalSet.end();

		while(signalIter != signalIterEnd)
		{
			const SignalPtr signal = *signalIter;
			jointDimension += signal->dimension();

			++signalIter;
		}

		const Integer2 sharedTime = 
			sharedTimeInterval(signalSet, lagSet);

		if (tMin)
		{
			*tMin = sharedTime[0];
		}

		const integer samples = sharedTime[1] - sharedTime[0];

		// Allocate the joint signal.

		SignalPtr jointSignal(new Signal(samples, jointDimension));
		
		if  (samples == 0)
		{
			// There is no common time interval that
			// all signals would share.
			return jointSignal;
		}

		// Copy the signals into parts of the joint signal.

		const integer maxLag = *std::max_element(
			lagSet.begin(), lagSet.end());

		integer dimensionOffset = 0;

		Integer_Iterator lagIter = lagSet.begin();
		const Integer_Iterator lagIterEnd = lagSet.end();

		signalIter = signalSet.begin();
		while(lagIter != lagIterEnd)
		{
			const SignalPtr signal = *signalIter;
			const integer lagOffset = maxLag - *lagIter;
			const integer dimension = signal->dimension();

			for (integer i = 0;i < samples;++i)
			{
				std::copy(
					signal->data().rowBegin(i + lagOffset),
					signal->data().rowEnd(i + lagOffset),
					jointSignal->data().rowBegin(i) + dimensionOffset);
			}

			dimensionOffset += dimension;

			++signalIter;
			++lagIter;
		}

		return jointSignal;
	}

	template <typename SignalPtr_Iterator>
	SignalPtr merge(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		integer* tMin)
	{
		return Tim::merge(signalSet,
			constantRange(0, signalSet.size()),
			tMin);
	}

	template <
		typename SignalPtr_OutputIterator,
		typename Integer_Iterator>
	void merge(
		const Array<SignalPtr, 2>& ensembleSet,
		SignalPtr_OutputIterator result,
		const ForwardRange<Integer_Iterator>& lagSet)
	{
		ENSURE_OP(lagSet.size(), ==, ensembleSet.height());

		const integer trials = ensembleSet.width();
		for (integer i = 0;i < trials;++i)
		{
			*result = merge(
				forwardRange(ensembleSet.columnBegin(i),
				ensembleSet.columnEnd(i)), lagSet);
			++result;
		}
	}

	template <typename SignalPtr_OutputIterator>
	void merge(
		const Array<SignalPtr, 2>& ensembleSet,
		SignalPtr_OutputIterator result)
	{
		Tim::merge(ensembleSet, result,
			constantRange(0, ensembleSet.height()));
	}

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_OutputIterator>
	void merge(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		SignalPtr_OutputIterator result,
		integer xLag,
		integer yLag)
	{
		ENSURE_OP(xSignalSet.size(), ==, ySignalSet.size());
		
		SignalPtr_X_Iterator xIter = xSignalSet.begin();
		const SignalPtr_X_Iterator xIterEnd = xSignalSet.end();
		SignalPtr_Y_Iterator yIter = ySignalSet.begin();

		while(xIter != xIterEnd)
		{
			*result = merge(*xIter, *yIter, xLag, yLag);
			
			++result;
			++xIter;
			++yIter;
		}
	}

}

#endif
