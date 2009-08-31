#include "mex.h"

#include "tim/core/transfer_entropy.h"

#include <pastel/sys/pastelomp.h>

#include <boost/static_assert.hpp>

#include <algorithm>

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
		xIndex,
		yIndex,
		zIndex,
		wIndex,
		yLagIndex,
		zLagIndex,
		wLagIndex,
		kNearestIndex,
		threadsIndex
	};

	const integer trials = mxGetNumberOfElements(inputSet[xIndex]);

	// X-signals

	std::vector<SignalPtr> xEnsemble;
	xEnsemble.reserve(trials);

	for (integer i = 0;i < trials;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[xIndex], i);

		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.
		const integer samples = mxGetN(signalArray);
		const integer dimension = mxGetM(signalArray);

		real* rawData = mxGetPr(signalArray);

		xEnsemble.push_back(SignalPtr(
			new Signal(samples, dimension, rawData)));
	}

	// Y-signals

	std::vector<SignalPtr> yEnsemble;
	yEnsemble.reserve(trials);

	for (integer i = 0;i < trials;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[yIndex], i);

		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.
		const integer samples = mxGetN(signalArray);
		const integer dimension = mxGetM(signalArray);

		real* rawData = mxGetPr(signalArray);

		yEnsemble.push_back(SignalPtr(
			new Signal(samples, dimension, rawData)));
	}

	// Z-signals

	std::vector<SignalPtr> zEnsemble;
	zEnsemble.reserve(trials);

	for (integer i = 0;i < trials;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[zIndex], i);

		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.
		const integer samples = mxGetN(signalArray);
		const integer dimension = mxGetM(signalArray);

		real* rawData = mxGetPr(signalArray);

		zEnsemble.push_back(SignalPtr(
			new Signal(samples, dimension, rawData)));
	}

	// W-signals

	std::vector<SignalPtr> wEnsemble;
	wEnsemble.reserve(trials);

	for (integer i = 0;i < trials;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[wIndex], i);

		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.
		const integer samples = mxGetN(signalArray);
		const integer dimension = mxGetM(signalArray);

		real* rawData = mxGetPr(signalArray);

		wEnsemble.push_back(SignalPtr(
			new Signal(samples, dimension, rawData)));
	}

	// Other parameters

	const integer yLag = *mxGetPr(inputSet[yLagIndex]);
	const integer zLag = *mxGetPr(inputSet[zLagIndex]);
	const integer wLag = *mxGetPr(inputSet[wLagIndex]);
	const integer kNearest = *mxGetPr(inputSet[kNearestIndex]);
	const integer threads = *mxGetPr(inputSet[threadsIndex]);

#if PASTEL_OMP != 0
	omp_set_num_threads(threads);
#endif

	outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	*rawResult = transferEntropy(
		forwardRange(xEnsemble.begin(), xEnsemble.end()),
		forwardRange(yEnsemble.begin(), yEnsemble.end()),
		forwardRange(zEnsemble.begin(), zEnsemble.end()),
		forwardRange(wEnsemble.begin(), wEnsemble.end()),
		yLag, zLag, wLag,
		kNearest);
}
