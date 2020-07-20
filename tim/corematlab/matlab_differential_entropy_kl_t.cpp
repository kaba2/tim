// Description: differential_entropy_kl_t
// DocumentationOf: differential_entropy_kl_t.m

#include "tim/corematlab/tim_matlab.h"

#include "tim/core/differential_entropy_kl.h"

void force_linking_differential_entropy_kl_t() {};

using namespace Tim;

namespace
{

	void matlabTemporalDifferentialEntropyKl(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			X,
			TimeWindowRadius,
			KNearest,
			FilterIndex,
			Inputs
		};

		enum Output
		{
			Estimate,
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		std::vector<MatlabMatrix<dreal>> xMatrices = matlabAsMatrixRange<dreal>(inputSet[X]) | ranges::to_vector;
		std::vector<Signal> xSignals = matlabMatricesAsSignals(xMatrices) | ranges::to_vector;

		integer timeWindowRadius = matlabAsScalar<integer>(inputSet[TimeWindowRadius]);
		integer kNearest = matlabAsScalar<integer>(inputSet[KNearest]);

		std::vector<dreal> filter;
		matlabGetScalars(inputSet[FilterIndex], std::back_inserter(filter));

		Signal estimate = temporalDifferentialEntropyKl(
			xSignals, 
			timeWindowRadius, 
			kNearest,
			Default_Norm(),
			filter);

		integer nans = std::max(estimate.t(), (integer)0);
		integer skip = std::max(-estimate.t(), (integer)0); 
		integer samples = std::max(nans + estimate.samples() - skip, (integer)0);

		MatrixView<dreal> result = matlabCreateMatrix<dreal>(1, samples, outputSet[Estimate]);
		ranges::fill(result.slicex(0, nans).range(), (dreal)Nan());
		ranges::copy(
			estimate.data().slicex(skip).range(),
			std::begin(result.slicex(nans).range()));
	}

	void addFunction()
	{
		matlabAddFunction(
			"differential_entropy_kl_t",
			matlabTemporalDifferentialEntropyKl);
	}

	CallFunction run(addFunction);

}
