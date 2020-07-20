// Description: Splitting a signal into lower-dimensional signals

#ifndef TIM_SIGNAL_SPLIT_H
#define TIM_SIGNAL_SPLIT_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

namespace Tim
{

	//! Creates aliases for 1d-marginal signals.
	/*!
	See the more specialized
	'split' function for more information.
	*/
	template <typename Signal_OutputIterator>
	void split(
		const Signal& jointSignal,
		Signal_OutputIterator signalSet)
	{
		integer dimension = jointSignal.dimension();

		std::vector<integer> partition;
		partition.reserve(dimension + 1);
		for (integer i = 0;i <= dimension;++i)
		{
			partition.push_back(i);
		}

		Tim::split(jointSignal, partition, signalSet);
	}

	//! Creates aliases for marginal signals.
	/*!
	Preconditions:
	partition[x] < partition[x + 1]

	The partition of the joint signal is given
	by the 'partition' set. Assume that 'jointSignal'
	is of dimension 4 and 'partition' contains the 
	numbers 0, 2, 3, 4. Then the marginal signals
	are split with the dimension subranges
	[0, 2[, [2, 3[, [3, 4[.
	*/
	template <typename Signal_OutputIterator>
	void split(
		const Signal& jointSignal,
		const std::vector<integer>& partition,
		Signal_OutputIterator signalSet)
	{
		ENSURE_OP(partition.size(), >=, 2);

		integer signals = partition.size() - 1;

		for (integer x = 0;x < signals;++x)
		{
			PENSURE_OP(partition[x], <, partition[x + 1]);

			integer marginalDimension = 
				partition[x + 1] - partition[x];


			*signalSet = Tim::split(jointSignal, partition[x], 
				partition[x] + marginalDimension);
			++signalSet;
		}
	}

	//! Creates an alias for a marginal signal.
	/*
	Returns a signal S which refers to a marginal signal of P. 
	For each sample x in P, S refers only to the dimension subrange 
	[dimensionBegin, dimensionEnd[. S then has dimension
	(dimensionEnd - dimensionBegin). Note S shares memory
	with P and thus changes in either are reflected in
	the other.
	*/
	inline TIM Signal split(
		const Signal& signal,
		integer dimensionBegin,
		integer dimensionEnd)
	{
		ENSURE_OP(dimensionBegin, <=, dimensionEnd);
		ENSURE_OP(dimensionBegin, >=, 0);
		ENSURE_OP(dimensionEnd, <=, signal.dimension());

		integer dimension = dimensionEnd - dimensionBegin;
		integer samples = signal.samples();

		return Signal(signal.data().slicex(dimensionBegin, dimensionEnd), signal.t());
	}

}

#endif
