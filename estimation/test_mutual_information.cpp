#include "estimation.h"

#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/mutual_information.h"

#include <pastel/math/matrix_tools.h>
#include <pastel/math/cholesky_decomposition_tools.h>

#include <pastel/sys/string_tools.h>
#include <pastel/sys/measuretable_tools.h>

#include <pastel/device/timer.h>

#include <iomanip>

using namespace Tim;

namespace
{

	void testMutualInformationCase(
		const std::string& name,
		const SignalPtr& jointSignal,
		integer kNearest,
		real maxRelativeError,
		real correctMi)
	{
		const real mi = mutualInformation(
			jointSignal,
			kNearest,
			maxRelativeError);

		log() << name << ": " << mi
			<< " (" << mi - correctMi << ")"
			<< logNewLine;
	}

	void testMutualInformation()
	{
		log() << "Mutual information estimates: " << logNewLine;

		const integer samples = 100000;
		const integer kNearest = 1;
		const real maxRelativeError = 0;

		log() << "2d correlated gaussian" << logNewLine;

		{
			const integer dimension = 2;

			for (integer i = 0;i < 10;++i)
			{
				MatrixD covariance(dimension, dimension);

				const real r = (real)i / 10;

				covariance |= 
					1, r,
					r, 1;

				const CholeskyDecompositionD cholesky(
					covariance);

				const real det = determinant(cholesky);
				const real cond = (1 + r) / (1 - r);

				ENSURE(cholesky.succeeded());

				const SignalPtr jointSignal = 
					generateCorrelatedGaussian(samples, dimension, cholesky);

				testMutualInformationCase(
					"Cor.Gauss. det " + realToString(det) +
					" cond " + realToString(cond),
					jointSignal,
					kNearest,
					maxRelativeError,
					correlatedGaussianMutualInformation(
					diagonalProduct(covariance), determinant(cholesky)));
			}
		}

		log() << "nD correlated gaussian cond-det covariances" << logNewLine;
		{
			const integer dimension = 5;
			for (integer i = 0;i < 10;++i)
			{
				MatrixD covariance(dimension, dimension);

				const real cond = 10 - i;
				const real det = 1 + i;

				//const real cond = 2;
				//const real det = 1 + i;

				/*
				const real cond = 1 + i;
				const real det = 2;
				*/

				/*
				setRandomMatrix(covariance);
				covariance *= transpose(covariance);
				*/

				setRandomSymmetricPositiveDefinite(
					det, cond, covariance);

				/*
				log() << "cond = " << conditionManhattan(covariance)
					<< ", det = " << determinant(covariance) << logNewLine;
				*/

				const CholeskyDecompositionD cholesky(
					covariance);

				ENSURE(cholesky.succeeded());

				std::cout << std::fixed << std::setprecision(4);
				
				/*
				std::cout << covariance << std::endl;

				std::cout << cholesky.lower() << std::endl;

				std::cout << determinant(cholesky) << ", " << det << std::endl;
				*/

				const SignalPtr jointSignal = 
					generateCorrelatedGaussian(samples, dimension, cholesky);

				//normalizeCovariance(jointSignal, covariance);

				/*
				MatrixD sampleCovariance;
				computeCovariance(jointSignal, sampleCovariance);
				std::cout << sampleCovariance << std::endl;
				*/

				testMutualInformationCase(
					"Cor.Gauss. det " + realToString(det) +
					" cond " + realToString(cond),
					jointSignal,
					kNearest,
					maxRelativeError,
					correlatedGaussianMutualInformation(
					diagonalProduct(covariance), determinant(cholesky)));
			}
		}
	}

	void testTiming()
	{
		log() << "Mutual information estimates: " << logNewLine;

		MeasureTable measureTable;
		measureTable.setCaption("Comparison between eeglab and tim.");
		measureTable.setSize(7, 21);

		enum
		{
			Samples_Column,
			Covariance_Column,
			ElTime_Column,
			TimTime_Column,
			ElMi_Column,
			TimMi_Column,
			CorrectMi_Column
		};

		measureTable(Samples_Column, 0).text() = "Samples";
		measureTable(Covariance_Column, 0).text() = "Cov.";
		measureTable(ElTime_Column, 0).text() = "El time";
		measureTable(TimTime_Column, 0).text() = "Tim time";
		measureTable(ElMi_Column, 0).text() = "El mi";
		measureTable(TimMi_Column, 0).text() = "Tim mi";
		measureTable(CorrectMi_Column, 0).text() = "Cor. mi";
		measureTable.addHorizontalSeparator(0);
		measureTable.addHorizontalSeparator(1);
		measureTable.addHorizontalSeparator(6);
		measureTable.addHorizontalSeparator(11);
		measureTable.addHorizontalSeparator(16);
		measureTable.addHorizontalSeparator(21);
		measureTable.addVerticalSeparator(Samples_Column);
		measureTable.addVerticalSeparator(ElTime_Column);
		measureTable.addVerticalSeparator(ElMi_Column);
		measureTable.addVerticalSeparator(CorrectMi_Column + 1);

		integer experiment = 1;
		const integer kNearest = 1;
		const integer dimension = 2;
		const real maxRelativeError = 0;
		for (integer i = 0;i < 4;++i)
		{
			const integer samples = 100 * (integer)std::pow((real)10, (real)i);

			for (integer j = 0;j < 5;++j)
			{
				MatrixD covariance(dimension, dimension);

				const real r = (real)j / 5;

				covariance |= 
					1, r,
					r, 1;

				const CholeskyDecompositionD cholesky(
					covariance);

				ENSURE(cholesky.succeeded());

				const SignalPtr jointSignal = 
					generateCorrelatedGaussian(samples, dimension, cholesky);

				const real correctMi = correlatedGaussianMutualInformation(
					diagonalProduct(covariance), determinant(cholesky));

				measureTable(Samples_Column, experiment).text() = 
					integerToString(samples);
				measureTable(Covariance_Column, experiment).text() =
					realToString(r, 4);
				measureTable(CorrectMi_Column, experiment).text() = 
					realToString(correctMi, 4);

				Timer timer;

				timer.setStart();
				const real mi = mutualInformation(
					jointSignal,
					kNearest,
					maxRelativeError);
				timer.store();

				measureTable(TimTime_Column, experiment).text() = 
					realToString(timer.seconds(), 4);
				measureTable(TimMi_Column, experiment).text() = 
					realToString(mi, 4);

				MatrixD pairwiseMi;
				timer.setStart();
				mutualInformationNaive(
					jointSignal,
					100,
					pairwiseMi);
				timer.store();

				measureTable(ElTime_Column, experiment).text() = 
					realToString(timer.seconds(), 4);
				measureTable(ElMi_Column, experiment).text() = 
					realToString(pairwiseMi(1, 0), 4);

				++experiment;
			}
		}

		printPretty(measureTable, std::cout);
		//printLatex(measureTable, std::cout);
	}

	void testAdd()
	{
		timTestList().add("mutual_information", testMutualInformation);
		timTestList().add("mi_timing", testTiming);
	}

	CallFunction run(testAdd);

}