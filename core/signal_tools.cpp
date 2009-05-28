#include "tim/core/signal_tools.h"

#include <pastel/math/cholesky_decomposition.h>
#include <pastel/math/uniformsampling.h>

#include <pastel/sys/view_all.h>

namespace Tim
{

	TIMCORE SignalPtr generateUniform(
		integer dimension,
		integer size)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(size >= 0, size);

		SignalPtr signal = SignalPtr(new Signal(dimension, size));

		for (integer i = 0;i < size;++i)
		{
			(*signal)[i] = asPoint(randomVector<Dynamic, real>(dimension));
		}

		return signal;
	}

	TIMCORE SignalPtr generateGaussian(
		integer dimension,
		integer size)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(size >= 0, size);

		SignalPtr signal = SignalPtr(new Signal(dimension, size));

		for (integer i = 0;i < size;++i)
		{
			(*signal)[i] = asPoint(randomVectorGaussian<Dynamic, real>(dimension));
		}

		return signal;
	}

	TIMCORE SignalPtr generateCorrelatedGaussian(
		integer dimension,
		integer size,
		const DynamicMatrix& correlation)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(size >= 0, size);

		SignalPtr correlatedGaussian = generateGaussian(dimension, size);

		const CholeskyDecomposition<Dynamic, real> cholesky(correlation);
		
		if (!cholesky.succeeded())
		{
			log() << "Correlation matrix not (numerically) positive definite! Can't generate correlated gaussians." 
				<< logNewLine;
			return correlatedGaussian;
		}
	
		for (integer i = 0;i < size;++i)
		{
			(*correlatedGaussian)[i] = asPoint((*correlatedGaussian)[i].asVector() * cholesky.lower());
		}

		return correlatedGaussian;
	}

	TIMCORE void splitDimensions(
		const SignalPtr& signal,
		std::vector<SignalPtr>& signalSet)
	{
		ENSURE(!signal.empty());

		const integer dimension = signal->dimension();
		const integer size = signal->size();

		std::vector<SignalPtr> result;
		result.reserve(dimension);
		for (integer i = 0;i < dimension;++i)
		{
			const SignalPtr smallAliasSignal = 
				SignalPtr(new Signal(signal, i, i + 1));
			result.push_back(smallAliasSignal);
		}
		
		result.swap(signalSet);
	}

	TIMCORE SignalPtr mergeDimensions(
		const std::vector<SignalPtr>& signalList)
	{
		if (signalList.empty())
		{
			return SignalPtr();
		}

		const integer size = signalList.front()->size();

		const integer signals = signalList.size();
		integer bigDimension = 0;
		for (integer i = 0;i < signals;++i)
		{
			bigDimension += signalList[i]->dimension();
			ENSURE2(signalList[0]->size() == signalList[i]->size(), 
				signalList[0]->size(), signalList[i]->size());
		}

		SignalPtr bigSignal(new Signal(bigDimension, size));
		
		integer dimensionOffset = 0;

		for (integer i = 0;i < signals;++i)
		{
			copy(signalList[i]->constView(),
				subView(bigSignal->view(), 
				Rectangle2(dimensionOffset, 0, 
				dimensionOffset + signalList[i]->dimension(), size)));

			dimensionOffset += signalList[i]->dimension();
		}

		return bigSignal;
	}

}
