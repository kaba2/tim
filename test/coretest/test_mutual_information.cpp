#include "estimation.h"

#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/mutual_information.h"

#include <pastel/math/matrix_algorithms.h>
#include <pastel/math/cholesky_decomposition.h>

#include <pastel/sys/string_algorithms.h>
#include <pastel/sys/measuretable_tools.h>

#include <iomanip>
#include <numeric>

using namespace Tim;

namespace
{

	void testMutualInformationCase(
		const std::string& name,
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		integer xLag,
		integer yLag,
		integer timeWindowRadius,
		integer kNearest,
		real correctMi)
	{
		const real mi = mutualInformation(
			constantRange(xSignal), constantRange(ySignal), 
			xLag, yLag, kNearest);

		const SignalPtr temporalMi = temporalMutualInformation(
			constantRange(xSignal), constantRange(ySignal), 
			timeWindowRadius,
			xLag, yLag, kNearest);
		const real averageMi = 
			std::accumulate(temporalMi->data().begin(), 
			temporalMi->data().end(), (real)0) / temporalMi->samples();

		/*
		std::vector<SignalPtr> signalSet;
		split(jointSignal, signalSet);

		const real mi = mutualInformationFromEntropy(
			signalSet,
			kNearest,
			Euclidean_NormBijection<real>());
		*/

		log() << name << logNewLine;
		log() << "  average " << mi << "(" << mi - correctMi << ")" << logNewLine;
		log() << "  temporal average " << averageMi << "(" << averageMi - correctMi << ")" << logNewLine;

		/*
		const real re = relativeError<real>(averageMi, correctMi);

		log() << name << ": " << averageMi
			<< " (de = " << averageMi - correctMi << ", " << re * 100 << "%)"
			<< logNewLine;
		*/
	}

	void testMutualInformation()
	{
		log() << "Mutual information estimates: " << logNewLine;

		const integer samples = 10000;
		const integer timeWindowRadius = 1000;
		const integer kNearest = 10;

		log() << "2d correlated gaussian" << logNewLine;

		{
			const integer dimension = 2;

			for (integer i = 0;i < 10;++i)
			{
				Matrix<real> covariance(dimension, dimension);

				const real r = (real)i / 10;

				covariance |= 
					1, r,
					r, 1;

				const CholeskyDecomposition<real> cholesky(
					covariance);

				const real det = determinant(cholesky);
				const real cond = (1 + r) / (1 - r);

				ENSURE(cholesky.succeeded());

				const SignalPtr jointSignal = 
					generateCorrelatedGaussian(samples, dimension, cholesky);

				SignalPtr xSignal = split(jointSignal, 0, 1);
				SignalPtr ySignal = split(jointSignal, 1, 2);

				testMutualInformationCase(
					"Cor.Gauss. det " + realToString(det) +
					" cond " + realToString(cond),
					xSignal,
					ySignal,
					0,
					0,
					timeWindowRadius,
					kNearest,
					correlatedGaussianMutualInformation(
					diagonalProduct(covariance), determinant(cholesky)));
			}
		}

		log() << "nD correlated gaussian cond-det covariances" << logNewLine;
		{
			const integer dimension = 2;
			for (integer i = 0;i < 10;++i)
			{
				Matrix<real> covariance(dimension, dimension);

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

				const CholeskyDecomposition<real> cholesky(
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
				const SignalPtr xSignal = split(jointSignal, 0, 1);
				const SignalPtr ySignal = split(jointSignal, 1, 2);

				//normalizeCovariance(jointSignal, covariance);

				/*
				Matrix<real> sampleCovariance;
				computeCovariance(jointSignal, sampleCovariance);
				std::cout << sampleCovariance << std::endl;
				*/

				testMutualInformationCase(
					"Cor.Gauss. det " + realToString(det) +
					" cond " + realToString(cond),
					xSignal,
					ySignal,
					0,
					0,
					timeWindowRadius,
					kNearest,
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
		for (integer i = 0;i < 4;++i)
		{
			const integer samples = 100 * (integer)std::pow((real)10, (real)i);

			for (integer j = 0;j < 5;++j)
			{
				Matrix<real> covariance(dimension, dimension);

				const real r = (real)j / 5;

				covariance |= 
					1, r,
					r, 1;

				const CholeskyDecomposition<real> cholesky(
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

				SignalPtr xSignal = split(jointSignal, 0, 1);
				SignalPtr ySignal = split(jointSignal, 1, 2);

				//Timer timer;
				
				//timer.setStart();
				const real mi = mutualInformation(
					constantRange(xSignal),
					constantRange(ySignal),
					0, 0,
					kNearest);
				//timer.store();

				measureTable(TimTime_Column, experiment).text() = 
					realToString(/*timer.seconds()*/0, 4);
				measureTable(TimMi_Column, experiment).text() = 
					realToString(mi, 4);
				
				//timer.setStart();
				Array<real> pairwiseMi = mutualInformationFromBinning(
					jointSignal,
					100);
				//timer.store();

				measureTable(ElTime_Column, experiment).text() = 
					realToString(/*timer.seconds()*/0, 4);
				measureTable(ElMi_Column, experiment).text() = 
					realToString(pairwiseMi(1, 0), 4);

				++experiment;
			}
		}

		printPretty(measureTable, std::cout);
		//printLatex(measureTable, std::cout);
	}

	void testBoundaryLag()
	{
		SignalPtr signal = generateGaussian(100, 1);
		const SignalPtr temporalMi = temporalMutualInformation(
			constantRange(signal), 
			constantRange(signal), 
			10, 
			0, 99, 1);
	}

	void testTemporal()
	{
		//SignalPtr signal(new Signal(10, 1));
		//copy_n(countingIterator(0), 10, signal->data().begin());
		SignalPtr signal = generateGaussian(10, 1);
		real filter[] = {0.25, 0.5, 1, 0.5, 0.25};
		const SignalPtr temporalMi = temporalMutualInformation(
			constantRange(signal, 2),
			constantRange(signal, 2),
			2,
			0, 0, 1,
			range(filter));
	}

	void testAdd()
	{
		timTestList().add("mutual_information", testMutualInformation);
		timTestList().add("mi_timing", testTiming);
		timTestList().add("mi_boundary_lag", testBoundaryLag);
		timTestList().add("mi_temporal", testTemporal);
	}

	CallFunction run(testAdd);

}