#include "tim/core/signal_tools.h"

#include <pastel/math/cholesky_decomposition.h>

#include <pastel/sys/random_vector.h>
#include <pastel/sys/view_all.h>

namespace Tim
{

	TIMCORE SignalPtr generateUniform(
		integer size,
		integer dimension)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(size >= 0, size);

		SignalPtr signal = SignalPtr(new Signal(size, dimension));

		for (integer i = 0;i < size;++i)
		{
			(*signal)[i] = asPoint(randomVectorCube<Dynamic, real>(dimension));
		}

		return signal;
	}

	TIMCORE SignalPtr generateGaussian(
		integer size,
		integer dimension)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(size >= 0, size);

		SignalPtr signal = SignalPtr(new Signal(size, dimension));

		for (integer i = 0;i < size;++i)
		{
			(*signal)[i] = asPoint(randomGaussianVector<Dynamic, real>(dimension));
		}

		return signal;
	}

	TIMCORE SignalPtr generateCorrelatedGaussian(
		integer size,
		integer dimension,
		const CholeskyDecomposition<Dynamic, real>& correlationCholesky)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(size >= 0, size);
		ENSURE(correlationCholesky.lower().width() == dimension);
		ENSURE(correlationCholesky.succeeded());

		SignalPtr correlatedGaussian = generateGaussian(size, dimension);
	
		for (integer i = 0;i < size;++i)
		{
			(*correlatedGaussian)[i] = asPoint((*correlatedGaussian)[i].asVector() * correlationCholesky.lower());
		}

		return correlatedGaussian;
	}

	TIMCORE SignalPtr generateGeneralizedGaussian(
		integer size,
		integer dimension,
		real shape,
		real scale)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(size >= 0, size);

		SignalPtr signal = SignalPtr(new Signal(size, dimension));

		for (integer i = 0;i < size;++i)
		{
			(*signal)[i] = asPoint(randomGeneralizedGaussianVector<Dynamic, real>(dimension, shape, scale));
		}

		return signal;
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
				SignalPtr(new Signal(signal, i, 1));
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

		SignalPtr bigSignal(new Signal(size, bigDimension));
		
		integer dimensionOffset = 0;

		for (integer i = 0;i < signals;++i)
		{
			copy(signalList[i]->constView(),
				subView(bigSignal->view(), 
				Rectangle2(0, dimensionOffset, size, 
				dimensionOffset + signalList[i]->dimension())));

			dimensionOffset += signalList[i]->dimension();
		}

		return bigSignal;
	}

}
