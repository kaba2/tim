#include "mex.h"

#include "tim/core/entropy_combination.h"

#include <boost/static_assert.hpp>

#include <pastel/sys/pastelomp.h>
#include <pastel/sys/stdext_copy_n.h>

using namespace Tim;

void mexFunction(int outputs, mxArray *outputSet[],
				 int inputs, const mxArray *inputSet[])
{
	enum
	{
		RealIsDouble = boost::is_same<real, double>::value
	};
	BOOST_STATIC_ASSERT(RealIsDouble);

	enum
	{
		signalSetIndex,
		rangeSetIndex,
		lagSetIndex,
		kNearestIndex,
		threadsIndex
	};

	// In the following, keep in mind that 
	// Matlab uses column-major storage
	// while we use row-major storage. That is
	// why some values are assigned "the wrong way".

	const mxArray* signalSetArray = inputSet[signalSetIndex];

	const integer signals = mxGetM(signalSetArray);
	const integer trials = mxGetN(signalSetArray);

	Array<SignalPtr, 2> signalSet(trials, signals);

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

	std::vector<integer> lagSet;
	lagSet.reserve(signals);
	StdExt::copy_n(mxGetPr(inputSet[lagSetIndex]), signals, std::back_inserter(lagSet));

	const integer marginals = mxGetM(inputSet[rangeSetIndex]);
	std::vector<Integer3> rangeSet;
	rangeSet.reserve(marginals);
	{
		real* rawData = mxGetPr(inputSet[rangeSetIndex]);
		for (integer i = 0;i < marginals;++i)
		{
			rangeSet.push_back(Integer3(*rawData - 1, *(rawData + marginals), *(rawData + 2 * marginals)));
			//printf("%d %d %d ", rangeSet.back()[0], rangeSet.back()[1], rangeSet.back()[1]);
			++rawData;
		}
	}

	const integer kNearest = *mxGetPr(inputSet[kNearestIndex]);
	const integer threads = *mxGetPr(inputSet[threadsIndex]);

#if PASTEL_ENABLE_OMP != 0
	omp_set_num_threads(threads);
#endif

	outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	*rawResult = entropyCombination(
		signalSet,
		forwardRange(rangeSet.begin(), rangeSet.end()),
		forwardRange(lagSet.begin(), lagSet.end()),
		kNearest);
}
