#ifndef TIM_DIFFERENTIAL_ENTROPY_SP_HPP
#define TIM_DIFFERENTIAL_ENTROPY_SP_HPP

#include "tim/core/differential_entropy_sp.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/vector/vector_tools.h>
#include <pastel/sys/math_functions.h>
#include <pastel/sys/sequence/copy_n.h>

#include <algorithm>

namespace Tim
{

	namespace Detail_DifferentialEntropySp
	{

		class Computation
		{
		public:
			Computation(
				integer samples,
				integer dimension,
				integer minLevel)
				: samples_(samples)
				, dimension_(dimension)
				, minLevel_(minLevel)
			{
				unused(dimension_);
			}

			template <typename Iterator>
			real work(
				const Iterator& begin,
				const Iterator& end,
				const VectorD& min,
				const VectorD& max,
				integer level) const
			{
				// "Fast Multidimensional Entropy Estimation
				// by k-d Partitioning", 
				// Dan Stowell, Mark D. Plumbley,
				// IEEE Signal Processing Letters, Vol. 16, No. 6,
				// June 2009.

				if (begin == end)
				{
					return 0;
				}

				integer n = end - begin;
				Iterator medianIter = begin + n / 2;
				integer d = maxIndex(max - min);

				Compare compare(d);
				std::nth_element(begin, medianIter, end, compare);


				const real median = (*medianIter)[d];
				
				const real z = 2 * std::sqrt((real)n) * 
					(median - linear(min[d], max[d], 0.5)) /
					(max[d] - min[d]);
				
				if (n <= 3 || (level >= minLevel_ && std::abs(z) >= 1.96))
				{
					real p = (real)n / samples_;


					return p * std::log(product(max - min) / p);
				}
				else
				{
					VectorD leftMax = max;
					leftMax[d] = median;

					VectorD rightMin = min;
					rightMin[d] = median;

					return work(begin, medianIter, min, leftMax, level + 1) +
						work(medianIter, end, rightMin, max, level + 1);
				}
			}
		
		private:
			class Compare
			{
			public:
				explicit Compare(integer dimension)
					: dimension_(dimension)
				{
				}

				bool operator()(
					const real* left,
					const real* right) const
				{
					if (left[dimension_] < right[dimension_])
					{
						return true;
					}

					if (right[dimension_] < left[dimension_])
					{
						return false;
					}
					
					return left < right;
				}

			private:
				integer dimension_;
			};

			integer samples_;
			integer dimension_;
			integer minLevel_;
		};

	}

	template <typename SignalPtr_Range>
	real differentialEntropySp(
		const SignalPtr_Range& signalSet)
	{
		if (signalSet.empty())
		{
			return 0;
		}

		integer signals = signalSet.size();
		integer samples = minSamples(signalSet);

		const integer n = samples * signals;

		// Gather the point set.

		std::vector<const real*> pointSet;
		pointSet.reserve(n);
		auto iter = signalSet.begin();
		auto iterEnd = signalSet.end();
		while(iter != iterEnd)
		{
			const Signal& signal = **iter;
			copy_n(
				signal.pointBegin(), samples,
				std::back_inserter(pointSet));
			
			++iter;
		}
		
		// Compute bounds.

		integer dimension = signalSet.front()->dimension();
		VectorD min(ofDimension(dimension), infinity<real>());
		VectorD max(ofDimension(dimension), -infinity<real>());
		for (integer i = 0;i < n;++i)
		{

			const real* point = pointSet[i];
			for (integer d = 0;d < dimension;++d)
			{
				real position = point[d];
				if (position < min[d])
				{
					min[d] = position;
				}
				if (position > max[d])
				{
					max[d] = position;
				}
			}
		}

		// Compute differential entropy.

		integer minLevel = std::ceil(log2<real>(n) / 2);
		Detail_DifferentialEntropySp::Computation
			computation(n, dimension, minLevel);

		return computation.work(
			pointSet.begin(), pointSet.end(),
			min, max, 0);
	}

}

#endif
