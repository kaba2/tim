#ifndef TIM_DIVERGENCE_WKV_HPP
#define TIM_DIVERGENCE_WKV_HPP

#include "tim/core/divergence_wkv.h"
#include "tim/core/signalpointset.h"

#include <pastel/geometry/search_nearest.h>
#include <pastel/geometry/nearestset/kdtree_nearestset.h>

#include <pastel/sys/indicator/predicate_indicator.h>

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

namespace Tim
{

	template <
		typename X_Signal_Range,
		typename Y_Signal_Range>
	dreal divergenceWkv(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet)
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

		integer xDimension = std::begin(xSignalSet)->dimension();
		integer yDimension = std::begin(ySignalSet)->dimension();

		ENSURE_OP(xDimension, ==, yDimension);

		typedef SignalPointSet::Point_ConstIterator Point_ConstIterator;

		// Construct point-sets.

		SignalPointSet xPointSet(xSignalSet);
		SignalPointSet yPointSet(ySignalSet);

		integer xSamples = xPointSet.samples();
		integer ySamples = yPointSet.samples();

		using Block = tbb::blocked_range<integer>;
		using Pair = std::pair<dreal, integer>;
		
		auto compute = [&](
			const Block& block,
			const Pair& start)
		{
			dreal estimate = start.first;
			integer acceptedSamples = start.second;
			for (integer i = block.begin(); i < block.end(); ++i)
			{
				// Find out the nearest neighbor in X for a point in X.

				Point_ConstIterator query =
					*(xPointSet.begin() + i);

				Vector<dreal> queryPoint(
					ofDimension(xDimension),
					withAliasing((dreal*)(query->point())));

				dreal xxDistance2 =
					(dreal)searchNearest(
						kdTreeNearestSet(xPointSet.kdTree()), 
						queryPoint,
						PASTEL_TAG(accept), predicateIndicator(query, NotEqualTo())
					).first;

				if (xxDistance2 > 0 && xxDistance2 < infinity<dreal>())
				{
					// Find out the nearest neighbor in Y for a point in X.

					dreal xyDistance2 =
						(dreal)searchNearest(kdTreeNearestSet(yPointSet.kdTree()), 
							queryPoint).first;
					
					if (xyDistance2 > 0 && xyDistance2 < infinity<dreal>())
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

		dreal estimate = 0;
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
			estimate *= (dreal)xDimension / (2 * acceptedSamples);
			estimate += std::log((dreal)ySamples / (xSamples - 1));
		}
		else
		{
			estimate = (dreal)Nan();
		}

		return estimate;
	}

}

#endif
