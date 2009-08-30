#include "mex.h"

#include "tim/core/mutual_information_ec.h"

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

	const integer signals = mxGetDimensions(inputSet[0])[0];
	const integer trials = mxGetDimensions(inputSet[0])[1];

	std::vector<SignalPtr> xEnsemble;
	xEnsemble.reserve(trials);

	for (integer i = 0;i < trials;++i)
	{
		mxArray* signalArray = mxGetCell(inputSet[0], i);

		// It is intentional to assign the width
		// and height the wrong way. The reason
		// is that Matlab uses column-major storage
		// while we use row-major storage.
		const mwSize samples = mxGetN(signalArray);
		const mwSize dimension = mxGetM(signalArray);

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
		const mwSize samples = mxGetN(signalArray);
		const mwSize dimension = mxGetM(signalArray);

		real* rawData = mxGetPr(signalArray);

		yEnsemble.push_back(SignalPtr(
			new Signal(samples, dimension, rawData)));
	}

	const integer timeWindowRadius = *mxGetPr(inputSet[2]);
	const integer yLag = *mxGetPr(inputSet[3]);
	const integer kNearest = *mxGetPr(inputSet[4]);
	const integer threads = *mxGetPr(inputSet[5]);

#if PASTEL_OMP != 0
	omp_set_num_threads(threads);
#endif

	std::vector<real> result;
	temporalMutualInformation(
		forwardRange(xEnsemble.begin(), xEnsemble.end()),
		forwardRange(yEnsemble.begin(), yEnsemble.end()),
		timeWindowRadius,
		std::back_inserter(result),
		yLag,
		kNearest);

	outputSet[0] = mxCreateDoubleMatrix(1, result.size(), mxREAL);
	real* rawResult = mxGetPr(outputSet[0]);

	std::copy(result.begin(), result.end(), rawResult);
}
