// Description: Mutual information between marginals of a normal distribution
// DocumentationOf: mutual_information_normal.m

#include "tim/matlab/tim_matlab.h"

#include "tim/core/mutual_information_normal.h"

void force_linking_mutual_information_normal() {};

using namespace Tim;

namespace
{

	void matlabMutualInformationNormal(
		int outputs, mxArray *outputSet[],
		int inputs, const mxArray *inputSet[])
	{
		enum Input
		{
			ProductOfCovarianceDeterminants,
			JointCovarianceDeterminant,
			Inputs
		};

		enum Output
		{
			MutualInformation,
			Outputs
		};

		ENSURE_OP(inputs, ==, Inputs);
		ENSURE_OP(outputs, ==, Outputs);

		real productDetCov =
			matlabAsScalar<real>(inputSet[ProductOfCovarianceDeterminants]);
		real jointDetCov = 
			matlabAsScalar<real>(inputSet[JointCovarianceDeterminant]);

		real* mutualInformation = 
			matlabCreateScalar<real>(outputSet[MutualInformation]);
		*mutualInformation = 
			mutualInformationNormal<real>(productDetCov, jointDetCov);
	}

	void addFunction()
	{
		matlabAddFunction(
			"mutual_information_normal",
			matlabMutualInformationNormal);
	}

	CallFunction run(addFunction);

}
