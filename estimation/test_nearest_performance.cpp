#include <boost/random.hpp>

#include <pastel/device/devicesystem.h>
#include <pastel/device/timer.h>

#include <pastel/geometry/pointkdtree_tools.h>
#include <pastel/geometry/all_nearest_neighbors.h>
#include <pastel/geometry/point_patterns.h>

#include <pastel/sys/tuplebase.h>
#include <pastel/sys/randomdistribution.h>

#include <pastel/gfx/pcx.h>
#include "pastel/gfx/drawing.h"
#include "pastel/gfx/imagegfxrenderer.h"
#include "pastel/gfx/gfxrenderer_tools.h"
#include "pastel/gfx/color_tools.h"

#include <pastel/sys/vector_tools.h>
#include <pastel/sys/log_all.h>
#include <pastel/sys/random.h>
#include <pastel/sys/string_tools.h>
#include <pastel/sys/config_tools.h>
#include <pastel/sys/configfile.h>

#include <pastel/math/uniform_sampling.h>

#include <ANN/ANN.h>

#include <cmath>

#include "estimation.h"

using namespace Pastel;
using namespace Tim;

namespace
{

	bool configIncludePastelKdTree = true;
	bool configIncludeAnnBdTree = true;
	bool configIncludeAnnKdTree = true;
	bool configIncludeBruteForce = true;
	real configMaxRelativeError = 0;
	bool configUseDynamicVectors = true;
	bool configUseMultiThreading = false;
	bool configOutputLatexTable = false;
	std::vector<integer> configPointsSet;
	integer configFixedPoints = 100;
	std::vector<integer> configNeighborsSet;
	integer configFixedNeighbors = 1;
	std::vector<integer> configDimensionSet;
	integer configFixedDimension = 2;

	std::string space = "\t";
	std::string endOfLine = "";

	Log& timLog()
	{
		static Log theTimLog;

		return theTimLog;
	}

	template <int N, typename Real, typename NormBijection>
	void checkConsistency(const std::string& name,
						  const std::vector<Point<N, Real> >& pointSet,
						  const Array<2, integer>& correctNearest,
						  const Array<2, integer>& approxNearest,
						  const PASTEL_NO_DEDUCTION(Real)& maxRelativeError,
						  const NormBijection& normBijection)
	{
		ENSURE(correctNearest.extent() == approxNearest.extent());

		const integer kNearest = correctNearest.width();
		const integer points = correctNearest.height();

		Real maxError = 0;
		integer differ = 0;
		integer incorrect = 0;
		integer aNotFound = 0;
		integer bNotFound = 0;
		for (integer y = 0;y < points;++y)
		{
			for (integer x = 0;x < kNearest;++x)
			{
				if (correctNearest(x, y) < 0)
				{
					++aNotFound;
				}
				if (approxNearest(x, y) < 0)
				{
					++bNotFound;
				}
				if (correctNearest(x, y) < 0 ||
					approxNearest(x, y) < 0)
				{
					continue;
				}

				if (approxNearest(x, y) != correctNearest(x, y))
				{
					++differ;
				}

				const Real correctDistance = normBijection.toNorm(
					distance2(pointSet[correctNearest(x, y)], pointSet[y], normBijection));
				const Real approxDistance = normBijection.toNorm(
					distance2(pointSet[approxNearest(x, y)], pointSet[y], normBijection));

				const Real currentError = mabs(correctDistance - approxDistance) / correctDistance;

				if (currentError > maxError)
				{
					maxError = currentError;
				}
				if (currentError > maxRelativeError)
				{
					++incorrect;
				}
			}
		}

		const real differPercent = 
			((real)differ * 100) / (points * kNearest);

		const real incorrectPercent = 
			((real)incorrect * 100) / (points * kNearest);

		if (incorrect > 0)
		{
			log() << incorrect << " incorrect results detected!" 
				<< logNewLine;
		}
		/*
		log() << logNewLine;
		log() << name << " consistency:" << logNewLine;
		log() << differPercent << "% differs." << logNewLine;
		log() << incorrectPercent << "% incorrect." << logNewLine;
		log() << maxError << " max relative error." << logNewLine;
		*/
		/*
		timLog() << aNotFound << " not found with a, "
			<< bNotFound << " not found with b." << logNewLine;
		*/
	}

	template <int N, typename Real, typename NormBijection>
	void allNearestNeighborsKdTreeAnn(
		const std::vector<Point<N, Real> >& pointSet,
		integer kNearest,
		const PASTEL_NO_DEDUCTION(Real)& maxDistance,
		const PASTEL_NO_DEDUCTION(Real)& maxRelativeError,
		const NormBijection& normBijection,
		Array<2, integer>& nearestArray)
	{
		ENSURE1(kNearest > 0, kNearest);
		ENSURE2(kNearest < pointSet.size(), kNearest, pointSet.size());
		ENSURE1(maxDistance > 0, maxDistance);
		ENSURE1(maxRelativeError >= 0, maxRelativeError);

		const integer points = pointSet.size();
		if (points == 0)
		{
			return;
		}

		const integer dimension = pointSet.front().size();
		
		ANNpointArray annPointSet = new ANNpoint[points];
		for (integer i = 0;i < points;++i)
		{
			annPointSet[i] = (real*)&pointSet[i][0];
		}
		
		ANNidxArray annNearestSet = new ANNidx[kNearest + 1];
		ANNdistArray annDistanceSet = new ANNdist[kNearest + 1];

		ANNkd_tree* kdTree = 
			new ANNkd_tree(
			annPointSet,		
			points,			
			dimension,
			16);

		for (integer i = 0;i < points;++i)
		{
			kdTree->annkSearch(
				annPointSet[i],
				kNearest + 1,
				annNearestSet,
				annDistanceSet,
				maxRelativeError);

			integer nearestIndex = 0;
			for (integer j = 0;j < kNearest + 1 && nearestIndex < kNearest;++j)
			{
				if (annNearestSet[j] != i)
				{
					nearestArray(nearestIndex, i) = annNearestSet[j];
					++nearestIndex;
				}
			}
		}

		delete kdTree;
		delete [] annPointSet;
		delete [] annNearestSet;
		delete [] annDistanceSet;
	}

	template <int N, typename Real, typename NormBijection>
	void allNearestNeighborsBdTreeAnn(
		const std::vector<Point<N, Real> >& pointSet,
		integer kNearest,
		const PASTEL_NO_DEDUCTION(Real)& maxDistance,
		const PASTEL_NO_DEDUCTION(Real)& maxRelativeError,
		const NormBijection& normBijection,
		Array<2, integer>& nearestArray)
	{
		ENSURE1(kNearest > 0, kNearest);
		ENSURE2(kNearest < pointSet.size(), kNearest, pointSet.size());
		ENSURE1(maxDistance > 0, maxDistance);
		ENSURE1(maxRelativeError >= 0, maxRelativeError);

		const integer points = pointSet.size();
		if (points == 0)
		{
			return;
		}

		const integer dimension = pointSet.front().size();
		
		ANNidxArray annNearestSet = new ANNidx[kNearest + 1];
		ANNdistArray annDistanceSet = new ANNdist[kNearest + 1];

		ANNpointArray annPointSet = new ANNpoint[points];
		for (integer i = 0;i < points;++i)
		{
			annPointSet[i] = (real*)&pointSet[i][0];
		}

		ANNbd_tree* bdTree = 
			new ANNbd_tree(
			annPointSet,		
			points,			
			dimension,
			16);

		for (integer i = 0;i < points;++i)
		{
			bdTree->annkSearch(
				annPointSet[i],
				kNearest + 1,
				annNearestSet,
				annDistanceSet,
				maxRelativeError);

			integer nearestIndex = 0;
			for (integer j = 0;j < kNearest + 1 && nearestIndex < kNearest;++j)
			{
				if (annNearestSet[j] != i)
				{
					nearestArray(nearestIndex, i) = annNearestSet[j];
					++nearestIndex;
				}
			}
		}

		delete bdTree;
		delete [] annPointSet;
		delete [] annNearestSet;
		delete [] annDistanceSet;
	}

	template <int N, typename Real>
	typename boost::disable_if_c<(N == 2)>::type
		drawNearest(const std::string& name,
		const std::vector<Point<N, Real> >& pointSet,
		const Array<2, integer>& neighborSet)
	{
	}

	template <int N, typename Real>
	typename boost::enable_if_c<(N == 2)>::type
		drawNearest(const std::string& name,
		const std::vector<Point<N, Real> >& pointSet,
		const Array<2, integer>& neighborSet)
	{
		if (pointSet.size() > 100)
		{
			return;
		}

		const integer points = neighborSet.height();
		const integer kNearest = neighborSet.width();

		Array<2, Color> image(768, 768);
		ImageGfxRenderer<Color> renderer;
		renderer.setImage(&image);
		renderer.setColor(Color(1));
		renderer.setViewWindow(AlignedBox2(-1, -1, 1, 1));
		renderer.setFilled(false);

		for (integer i = 0;i < points;++i)
		{
			if (neighborSet(0, i) >= 0 && neighborSet(0, i) < pointSet.size())
			{
				renderer.setColor(Color(1));
			}
			else
			{
				renderer.setColor(Color(0, 1, 0));
			}
			drawCircle(renderer, Sphere2(pointSet[i], 0.02));
			for (integer j = 0;j < kNearest;++j)
			{
				const integer neighbor = neighborSet(j, i);
				if (neighbor >= 0)
				{
					renderer.setColor(hsvToRgb(Color((real)j / kNearest, 1, 1)));
					drawArrow(renderer, Segment2(pointSet[i], pointSet[neighbor]), 0.04 / (j + 1));
				}
			}
		}

		savePcx(image, "nearest_" + name + ".pcx");
	}

	real report(real timing, real referenceTiming)
	{
		//return (integer)(100 * (timing / referenceTiming) + 0.5);
		return timing;
	}

	template <int N, typename Real>
	void testIt(integer dimension, integer points, integer kNearest, const Real& maxRelativeError)
	{
		std::vector<Point<N, Real> > pointSet;
		//generateGaussianPointSet(points, dimension, pointSet);
		//generateUniformBallPointSet(points, dimension, pointSet);
		//generateClusteredPointSet(points, dimension, 10, pointSet);
		//generateUniformCubePointSet(points, dimension, pointSet);

		CountedPtr<Clustered_RandomDistribution<N, Real> >
			clusteredDistribution = clusteredRandomDistribution<N, Real>(dimension);

		const integer clusters = 10;
		for (integer i = 0;i < clusters;++i)
		{
			clusteredDistribution->add(
				translate(
				scale(
				gaussianRandomDistribution<N, Real>(dimension), 
				evaluate(randomVector<N, Real>(dimension) * 0.05)),
				randomVector<N, Real>(dimension)));
		}

		CountedPtr<RandomDistribution<N, Real> >
			randomDistribution = clusteredDistribution;

		/*
		CountedPtr<RandomDistribution<N, Real> >
			randomDistribution = 
			translate(
			scale(
			gaussianRandomDistribution<N, Real>(dimension), 
			randomVector<N, Real>(dimension)),
			randomVector<N, Real>(dimension));
		*/

		for (integer i = 0;i < points;++i)
		{
			pointSet.push_back(
				randomDistribution->sample());
		}

		// Pack the points for nice locality.

		/*
		Array<2, Real> pointArray(dimension, points);
		std::vector<Point<N, Real> > pointSet;
		pointSet.reserve(points);
		for (integer i = 0;i < points;++i)
		{

			pointSet.push_back(Point<N, Real>(ofDimension(0)));
	#if TIM_DYNAMIC != 0
			pointSet.back() = TemporaryPoint<N, Real>(
				ofDimension(dimension),
				withAliasing(&pointArray(0, i)));
	#endif
			pointSet.back() = gpPointSet[i];
			ENSURE(N != Dynamic || &(pointSet.back()[0]) == &pointArray(0, i));
		}
		*/

		EuclideanNormBijection<Real> normBijection;
		//InfinityNormBijection<Real> normBijection;
		//ManhattanNormBijection<Real> normBijection;
		//MinkowskiNormBijection<Real> normBijection(1.5);

		Array<2, integer> bruteNearest(kNearest, pointSet.size());

		timLog() << dimension << space << points << space << kNearest;

		// Brute force

		real bruteTime = 1;

		if (configIncludeBruteForce)
		{
			Timer timer;

			timer.setStart();

			allNearestNeighborsBruteForce(
				pointSet, kNearest, infinity<Real>(),
				normBijection,
				bruteNearest);

			timer.store();

			bruteTime = timer.seconds();
			drawNearest("brute", pointSet, bruteNearest);

			timLog() << space << bruteTime;
		}


		// Pastel kd-tree

		if (configIncludePastelKdTree)
		{
			Array<2, integer> kdNearest(kNearest, pointSet.size());
			{
				Timer timer;

				timer.setStart();

				allNearestNeighborsKdTree(
					pointSet,
					0,
					kNearest, 
					infinity<Real>(), 
					maxRelativeError,
					normBijection,
					16,
					SlidingMidpoint2_SplitRule(),
					&kdNearest);

				timer.store();

				timLog() << space << report(timer.seconds(), bruteTime);
			}

			if (configIncludeBruteForce)
			{
				checkConsistency("pkdtree" + integerToString(maxRelativeError), 
					pointSet, bruteNearest, kdNearest, maxRelativeError, normBijection);
			}
			drawNearest("pkdtree" + integerToString(maxRelativeError), pointSet, kdNearest);
		}

		if (configIncludeAnnKdTree)
		{
			// ANN kd-tree
			Array<2, integer> kdNearest(kNearest, pointSet.size());
			{
				Timer timer;

				timer.setStart();

				allNearestNeighborsKdTreeAnn(
					pointSet, kNearest, infinity<Real>(), maxRelativeError,
					normBijection,
					kdNearest);

				timer.store();

				timLog() << space << report(timer.seconds(), bruteTime);
			}

			if (configIncludeBruteForce)
			{
				checkConsistency("akdtree" + integerToString(maxRelativeError),
					pointSet, bruteNearest, kdNearest, maxRelativeError, normBijection);
			}
			drawNearest("akdtree" + integerToString(maxRelativeError), pointSet, kdNearest);
		}

		if (configIncludeAnnBdTree)
		{
			// ANN bd-tree
			Array<2, integer> bdNearest(kNearest, pointSet.size());
			{
				Timer timer;

				timer.setStart();

				allNearestNeighborsBdTreeAnn(
					pointSet, kNearest, infinity<Real>(), maxRelativeError,
					normBijection,
					bdNearest);

				timer.store();

				timLog() << space << report(timer.seconds(), bruteTime);
			}

			if (configIncludeBruteForce)
			{
				checkConsistency("bkdtree" + integerToString(maxRelativeError), 
					pointSet, bruteNearest, bdNearest, maxRelativeError, normBijection);
			}
			drawNearest("bkdtree" + integerToString(maxRelativeError), pointSet, bdNearest);
		}

		timLog() << endOfLine << logNewLine;
	}

	void test(integer dimension, integer points, integer kNearest, 
			  bool fixed, real maxRelativeError)
	{
		if (!fixed)
		{
			testIt<Dynamic, real>(dimension, points, kNearest, maxRelativeError);
		}
		else
		{
			switch(dimension)
			{
				case 1:
					testIt<1, real>(dimension, points, kNearest, maxRelativeError);
					break;
				case 2:
					testIt<2, real>(dimension, points, kNearest, maxRelativeError);
					break;
				case 4:
					testIt<4, real>(dimension, points, kNearest, maxRelativeError);
					break;
				case 8:
					testIt<8, real>(dimension, points, kNearest, maxRelativeError);
					break;
				case 16:
					testIt<16, real>(dimension, points, kNearest, maxRelativeError);
					break;
				case 32:
					testIt<32, real>(dimension, points, kNearest, maxRelativeError);
					break;
			};
		}
	}

	/*
	Dimension range
	2, 4, 8, 16, 32

	Point count range
	5000, 10000, 20000, 40000, 80000

	Neighbor count range
	1, 2, 4, 8, 16

	Relative error range
	0, 1, 2, 3, 4
	*/



	void timings()
	{
		timLog() 
			<< "Dim" << space
			<< "Points" << space
			<< "Nbors" << space;
		
		if (configIncludeBruteForce)
		{
			timLog() << "Brute" << space;
		}
		if (configIncludePastelKdTree)
		{
			timLog() << "PKd" << space;
		}
		if (configIncludeAnnKdTree)
		{
			timLog() << "AKd" << space;
		}
		if (configIncludeAnnBdTree)
		{
			timLog() << "ABd" << space;
		}

		timLog() << logNewLine;

		//test(8, 10000, 4, !configUseDynamicVectors, maxRelativeError);

		// Vary each dimension independently,
		// while fixing all other dimensions.

		for (integer i = 0;i < configDimensionSet.size();++i)
		{
			test(
				configDimensionSet[i],
				configFixedPoints,
				configFixedNeighbors,
				!configUseDynamicVectors, 
				configMaxRelativeError);
		}

		for (integer i = 0;i < configPointsSet.size();++i)
		{
			test(
				configFixedDimension,
				configPointsSet[i],
				configFixedNeighbors,
				!configUseDynamicVectors, 
				configMaxRelativeError);
		}
		for (integer i = 0;i < configNeighborsSet.size();++i)
		{
			test(
				configFixedDimension,
				configFixedPoints,
				configNeighborsSet[i], 
				!configUseDynamicVectors, 
				configMaxRelativeError);
		}
	}

	void estimation()
	{
		Config config;

		if (!loadConfig("config_performance.txt", config, LoadConfig_Echo::OnlyFinite))
		{
			log() << "Could not load config_performance.txt. Aborting..." 
				<< logNewLine;
			return;
		}

		const char *requiredNameList[] = 
		{
			"include_pastel_kdtree",
			"include_ann_kdtree",
			"include_ann_bdtree",
			"include_brute_force",
			"max_relative_error",
			"use_dynamic_vectors",
			"use_multi_threading",
			"output_latex_table",
			"points_set",
			"fixed_points",
			"neighbors_set",
			"fixed_neighbors",
			"dimension_set",
			"fixed_dimension"
		};

		const integer requiredNames = 
			sizeof(requiredNameList) / sizeof(const char*);

		std::vector<std::string> requiredSet;
		createPropertyList(requiredNameList, 
			requiredNames, requiredSet);

		if (!checkCreated(config, requiredSet, true))
		{
			return;
		}

		configIncludePastelKdTree = (config.property<integer>("include_pastel_kdtree") != 0);
		configIncludeAnnKdTree = (config.property<integer>("include_ann_kdtree") != 0);
		configIncludeAnnBdTree = (config.property<integer>("include_ann_bdtree") != 0);
		configIncludeBruteForce = (config.property<integer>("include_brute_force") != 0);
		configMaxRelativeError = config.property<real>("max_relative_error");
		configUseDynamicVectors = (config.property<integer>("use_dynamic_vectors") != 0);
		configUseMultiThreading = (config.property<integer>("use_multi_threading") != 0);
		configOutputLatexTable = (config.property<integer>("output_latex_table") != 0);
		configPointsSet = config.propertyList<integer>("points_set");
		configFixedPoints = config.property<integer>("fixed_points");
		configNeighborsSet = config.propertyList<integer>("neighbors_set");
		configFixedNeighbors = config.property<integer>("fixed_neighbors");
		configDimensionSet = config.propertyList<integer>("dimension_set");
		configFixedDimension = config.property<integer>("fixed_dimension");
		
		if (configOutputLatexTable)
		{
			space = " & ";
			endOfLine = " \\\\";
		}
		else
		{
			space = "\t";
			endOfLine = "";
		}

		timings();
	}

	void addTest()
	{
		timLog().addObserver(LogObserverPtr(new StreamLogObserver(&std::cout)));
		timLog().addObserver(LogObserverPtr(new FileLogObserver("timlog.txt")));

		timTestList().add("nearest_performance", estimation);
	}

	CallFunction run(addTest);

}

