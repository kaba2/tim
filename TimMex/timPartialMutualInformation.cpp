#include "mex.h"

#include "tim/core/partial_mutual_information.h"

#include <boost/static_assert.hpp>

#include <pastel/sys/pastelomp.h>

using namespace Tim;

void mexFunction(int outputs, mxArray *outputSet[],
				 int inputs, const mxArray *inputSet[])
{
	enum
	{
		RealIsDouble = boost::is_same<real, double>::value
	};
	BOOST_STATIC_ASSERT(RealIsDouble);

	const integer trials = mxGetNumberOfElements(inputSet[0]);

	std::vector<SignalPtr> xEnsemble;
	xEnsemble.reserve(trials);

	for (integer i = 0;i < trials;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[0], i);

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

	std::vector<SignalPtr> yEnsemble;
	yEnsemble.reserve(trials);

	for (integer i = 0;i < trials;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[1], i);

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

	std::vector<SignalPtr> zEnsemble;
	zEnsemble.reserve(trials);

	for (integer i = 0;i < trials;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[2], i);

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

	const integer yLag = *mxGetPr(inputSet[3]);
	const integer zLag = *mxGetPr(inputSet[4]);
	const integer kNearest = *mxGetPr(inputSet[5]);
	const integer threads = *mxGetPr(inputSet[6]);

#if PASTEL_OMP != 0
	omp_set_num_threads(threads);
#endif

	outputSet[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	*rawResult = partialMutualInformation(
		forwardRange(xEnsemble.begin(), xEnsemble.end()),
		forwardRange(yEnsemble.begin(), yEnsemble.end()),
		forwardRange(zEnsemble.begin(), zEnsemble.end()),
		yLag,
		zLag,
		kNearest);
}
