#ifndef TIM_DIVERGENCE_WKV_HPP
#define TIM_DIVERGENCE_WKV_HPP

#include "tim/core/divergence_wkv.h"
#include "tim/core/signalpointset.h"

#include <pastel/geometry/pointkdtree/pointkdtree_search_nearest.h>
#include <pastel/sys/indicator/predicate_indicator.h>
#include <pastel/sys/predicate/notequalto.h>

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

namespace Tim
{

	template <
		typename X_SignalPtr_Range,
		typename Y_SignalPtr_Range>
	real divergenceWkv(
		const X_SignalPtr_Range& xSignalSet,
		const Y_SignalPtr_Range& ySignalSet)
	{
		// "A Nearest-Neighbor Approach to Estimating
		// Divergence between Continuous Random Vectors"
		// Qing Wang, Sanjeev R. Kulkarni, Sergio Verdu,
		// IEEE International Symposium on Information Theory (ISIT), 
		// 2006.

		if (xSignalSet.empty() || ySignalSet.empty())
		{
			return 0;
		}

		integer xDimension = xSignalSet.front()->dimension();
		integer yDimension = ySignalSet.front()->dimension();

		ENSURE_OP(xDimension, ==, yDimension);

		typedef SignalPointSet::Point_ConstIterator Point_ConstIterator;

		// Construct point-sets.

		SignalPointSet xPointSet(xSignalSet);
		SignalPointSet yPointSet(ySignalSet);

		integer xSamples = xPointSet.samples();
		integer ySamples = yPointSet.samples();

		using Block = tbb::blocked_range<integer>;
		using Pair = std::pair<real, integer>;
		
		auto compute = [&](
			const Block& block,
			const Pair& start)
		{
			real estimate = start.first;
			integer acceptedSamples = start.second;
			for (integer i = block.begin(); i < block.end(); ++i)
			{
				// Find out the nearest neighbor in X for a point in X.

				Point_ConstIterator query =
					*(xPointSet.begin() + i);

				real xxDistance2 =
					searchNearest(xPointSet.kdTree(), query,
					nullOutput(),
					predicateIndicator(query, NotEqualTo()));


				if (xxDistance2 > 0 && xxDistance2 < infinity<real>())
				{
					// Find out the nearest neighbor in Y for a point in X.

					Vector<real> queryPoint(
						ofDimension(xDimension),
						withAliasing((real*)(query)->point()));

					real xyDistance2 =
						searchNearest(yPointSet.kdTree(), queryPoint);


					if (xyDistance2 > 0 && xyDistance2 < infinity<real>())
					{
						estimate += std::log(xyDistance2 / xxDistance2);
						++acceptedSamples;
					}
				}
			}

			return Pair(estimate, acceptedSamples);
		};

		auto reduce = [](const Pair& left, const Pair& right)
		{
			return Pair(
				left.first + right.first,
				left.second + right.second);
		};

		real estimate = 0;
		integer acceptedSamples = 0;

		std::tie(estimate, acceptedSamples) =
			tbb::parallel_reduce(
			Block(0, xSamples),
			Pair(0, 0),
			compute,
			reduce);

		if (acceptedSamples > 0)
		{
			// The factor 2 in the denominator is because 
			// 'xyDistance' and 'xxDistance' are squared distances
			// and thus need to be taken a square root. However,
			// this can be taken outside the logarithm with a 
			// division by 2.
			estimate *= (real)xDimension / (2 * acceptedSamples);
			estimate += std::log((real)ySamples / (xSamples - 1));
		}
		else
		{
			estimate = nan<real>();
		}

		return estimate;
	}

}

#endif
