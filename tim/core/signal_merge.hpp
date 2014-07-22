#ifndef TIM_SIGNAL_MERGE_HPP
#define TIM_SIGNAL_MERGE_HPP

#include "tim/core/signal_merge.h"
#include "tim/core/signal_properties.h"

namespace Tim
{

	template <
		typename Signal_Iterator,
		typename Integer_Iterator>
	Signal merge(
		const boost::iterator_range<Signal_Iterator>& signalSet,
		const boost::iterator_range<Integer_Iterator>& lagSet)
	{
		ENSURE_OP(signalSet.size(), ==, lagSet.size());

		if (signalSet.empty() ||
			lagSet.empty())
		{
			return Signal();
		}

		// Compute joint dimension.

		integer jointDimension = 0;
		Signal_Iterator signalIter = signalSet.begin();
		const Signal_Iterator signalIterEnd = signalSet.end();

		while(signalIter != signalIterEnd)
		{
			const Signal signal = *signalIter;
			jointDimension += signal.dimension();

			++signalIter;
		}

		const Integer2 sharedTime = 
			sharedTimeInterval(signalSet, lagSet);
		const integer samples = sharedTime[1] - sharedTime[0];

		// Allocate the joint signal.

		Signal jointSignal(samples, jointDimension, sharedTime[0]);
		
		if  (samples == 0)
		{
			// There is no common time interval that
			// all signals would share.
			return jointSignal;
		}

		// Copy the signals into parts of the joint signal.

		integer dimensionOffset = 0;

		Integer_Iterator lagIter = lagSet.begin();
		const Integer_Iterator lagIterEnd = lagSet.end();

		signalIter = signalSet.begin();
		while(lagIter != lagIterEnd)
		{
			const Signal signal = *signalIter;
			const integer lagOffset = sharedTime[0] - (signal.t() + *lagIter);
			const integer dimension = signal.dimension();

			for (integer i = 0;i < samples;++i)
			{
				std::copy(
					signal.data().cRowBegin(i + lagOffset),
					signal.data().cRowEnd(i + lagOffset),
					jointSignal.data().rowBegin(i) + dimensionOffset);
			}

			dimensionOffset += dimension;

			++signalIter;
			++lagIter;
		}

		return jointSignal;
	}

	template <typename Signal_Iterator>
	Signal merge(
		const boost::iterator_range<Signal_Iterator>& signalSet)
	{
		return Tim::merge(signalSet,
			constantRange(0, signalSet.size()));
	}

	template <
		typename Signal_OutputIterator,
		typename Integer_Iterator>
	void merge(
		const Array<Signal>& ensembleSet,
		Signal_OutputIterator result,
		const boost::iterator_range<Integer_Iterator>& lagSet)
	{
		ENSURE_OP(lagSet.size(), ==, ensembleSet.height());

		const integer trials = ensembleSet.width();
		for (integer i = 0;i < trials;++i)
		{
			*result = merge(
				range(ensembleSet.cColumnBegin(i),
				ensembleSet.cColumnEnd(i)), lagSet);
			++result;
		}
	}

	template <typename Signal_OutputIterator>
	void merge(
		const Array<Signal>& ensembleSet,
		Signal_OutputIterator result)
	{
		Tim::merge(ensembleSet, result,
			constantRange(0, ensembleSet.height()));
	}

	template <
		typename Signal_X_Iterator,
		typename Signal_Y_Iterator,
		typename Signal_OutputIterator>
	void merge(
		const boost::iterator_range<Signal_X_Iterator>& xSignalSet,
		const boost::iterator_range<Signal_Y_Iterator>& ySignalSet,
		Signal_OutputIterator result,
		integer xLag,
		integer yLag)
	{
		ENSURE_OP(xSignalSet.size(), ==, ySignalSet.size());
		
		Signal_X_Iterator xIter = xSignalSet.begin();
		const Signal_X_Iterator xIterEnd = xSignalSet.end();
		Signal_Y_Iterator yIter = ySignalSet.begin();

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
