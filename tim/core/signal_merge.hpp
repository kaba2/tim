#ifndef TIM_SIGNAL_MERGE_HPP
#define TIM_SIGNAL_MERGE_HPP

#include "tim/core/signal_merge.h"

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

		// In the following we think of having signals
		// embedded on the time axis, delayed with the given
		// lags. We accept to the merged signal only that
		// part in which all signals are present. E.g,
		// if the following are the time intervals that
		// three signals span:
		// 
		//        +--------------+
		//   +--------+
		//         +------+
		//
		// Then the their merged signal spans the following
		// time interval:
		//
		//         +--+

		// Find out the common time interval
		// [tLeftMax, tRightMin].

		integer samples = 0;
		integer maxLag = 0;
		{
			Integer_Iterator lagIter = lagSet.begin();
			const Integer_Iterator lagIterEnd = lagSet.end();
			SignalPtr_Iterator signalIter = signalSet.begin();

			integer tLeftMax = (*lagIter);
			integer tRightMin = tLeftMax + (*signalIter)->samples();

			while(lagIter != lagIterEnd)
			{
				const integer lag = *lagIter;

				if (lag > maxLag)
				{
					maxLag = lag;
				}

				const integer tLeft = lag;
				const integer tRight = lag + (*signalIter)->samples();

				if (tLeft > tLeftMax)
				{
					tLeftMax = tLeft;
				}
				if (tRight < tRightMin)
				{
					tRightMin = tRight;
				}

				++lagIter;
				++signalIter;
			}
			
			samples = tRightMin - tLeftMax;

			if (tMin)
			{
				*tMin = tLeftMax;
			}
		}

		if  (samples <= 0)
		{
			// There is no common time interval that
			// all signals would share.
			return SignalPtr(new Signal(0, jointDimension));
		}

		// Allocate the joint signal.

		SignalPtr jointSignal(new Signal(samples, jointDimension));
		
		// Copy the signals into parts of the joint signal.

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
			forwardRange(constantIterator(0)),
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
			forwardRange(constantIterator(0), ensembleSet.height()));
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
