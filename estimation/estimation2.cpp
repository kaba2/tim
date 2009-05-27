#define TIM_MULTI_THREADED 0
#define TIM_DYNAMIC 1
#define TIM_LATEX_OUTPUT 0
//#define TIM_ANN_BDTREE 1
//#define TIM_ANN_KDTREE 1
#define TIM_BRUTE_FORCE 1
#define TIM_KDTREE 1

#include <boost/random.hpp>

#include <pastel/device/devicesystem.h>
#include <pastel/device/timer.h>

#include <pastel/geometry/kdtree_tools.h>
#include <pastel/geometry/all_nearest_neighbors.h>

#include <pastel/sys/tuplebase.h>

#include <pastel/gfx/pcx.h>
#include "pastel/gfx/drawing.h"
#include "pastel/gfx/imagegfxrenderer.h"
#include "pastel/gfx/gfxrenderer_tools.h"
#include "pastel/gfx/color_tools.h"

#include <pastel/sys/vector_tools.h>
#include <pastel/sys/log_all.h>
#include <pastel/sys/random.h>
#include <pastel/sys/string_tools.h>

#include <pastel/math/uniformsampling.h>

#include <ANN/ANN.h>

#include <cmath>

using namespace Pastel;

#if TIM_LATEX_OUTPUT != 0
std::string space = " & ";
std::string endOfLine = " \\\\";
#else
std::string space = "\t";
std::string endOfLine = "";
#endif

Log& timLog()
{
	static Log theTimLog;

	return theTimLog;
}

template <int N, typename Real>
void generateUniformBallSet(integer points,
					 integer dimension,
					 std::vector<Point<N, Real> >& pointSet)
{
	std::vector<Point<N, Real> > result;
	result.reserve(points);

	for (integer i = 0;i < points;++i)
	{
		result.push_back(
			asPoint(randomVectorBall<N, Real>(dimension)));
	}

	result.swap(pointSet);
}

template <int N, typename Real>
void generateUniformCubeSet(integer points,
					 integer dimension,
					 std::vector<Point<N, Real> >& pointSet)
{
	std::vector<Point<N, Real> > result;
	result.reserve(points);

	for (integer i = 0;i < points;++i)
	{
		result.push_back(
			asPoint(randomVectorCube<N, Real>(dimension)));
	}

	result.swap(pointSet);
}

template <int N, typename Real>
void generateGaussianSet(integer points,
					 integer dimension,
					 std::vector<Point<N, Real> >& pointSet)
{
	std::vector<Point<N, Real> > result;
	result.reserve(points);

	for (integer i = 0;i < points;++i)
	{
		result.push_back(
			asPoint(randomVectorGaussian<N, Real>(dimension)));
	}

	result.swap(pointSet);
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

	log() << logNewLine;
	log() << name << " consistency:" << logNewLine;
	log() << differPercent << "% differs." << logNewLine;
	log() << incorrectPercent << "% incorrect." << logNewLine;
	log() << maxError << " max relative error." << logNewLine;
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
		dimension);

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
		4);

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
	generateGaussianSet(points, dimension, pointSet);
	//generateUniformBallSet(points, dimension, gpPointSet);
	//generateUniformCubeSet(points, dimension, gpPointSet);

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
#if TIM_BRUTE_FORCE != 0
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
	}

	timLog() << space << bruteTime;

#endif

	// Pastel kd-tree

#if TIM_KDTREE != 0
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
				&kdNearest);

			timer.store();

			timLog() << space << report(timer.seconds(), bruteTime);
		}

		checkConsistency("pkdtree" + integerToString(maxRelativeError), 
			pointSet, bruteNearest, kdNearest, maxRelativeError, normBijection);
		drawNearest("pkdtree" + integerToString(maxRelativeError), pointSet, kdNearest);
	}
#endif

#if TIM_ANN_KDTREE != 0
	// ANN kd-tree
	{
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

		checkConsistency("akdtree" + integerToString(maxRelativeError),
			pointSet, bruteNearest, kdNearest, maxRelativeError, normBijection);
		drawNearest("akdtree" + integerToString(maxRelativeError), pointSet, kdNearest);
	}
#endif

#if TIM_ANN_BDTREE != 0
	// ANN bd-tree
	{
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

		checkConsistency("bkdtree" + integerToString(maxRelativeError), 
			pointSet, bruteNearest, bdNearest, maxRelativeError, normBijection);
		drawNearest("bkdtree" + integerToString(maxRelativeError), pointSet, bdNearest);
	}
#endif

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
#if TIM_DYNAMIC == 0
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
#endif
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
	const bool multiThreaded = (TIM_MULTI_THREADED != 0 && PASTEL_ENABLE_OMP != 0);
	const bool fixed = (TIM_DYNAMIC == 0);
	const real maxRelativeError = 3;

	if (!multiThreaded)
	{
#if PASTEL_ENABLE_OMP != 0
		omp_set_num_threads(1);
#endif
		timLog() << "Single-threaded." << logNewLine;
	}
	else
	{
		timLog() << "Multi-threaded." << logNewLine;
	}

	if (fixed)
	{
		timLog() << "Static dimension." << logNewLine;
	}
	else
	{
		timLog() << "Dynamic dimension." << logNewLine;
	}

	timLog() << "Max relative error " << maxRelativeError 
		<< logNewLine;

	timLog() 
		<< "Dim" << space
		<< "Points" << space
		<< "Nbors" << space
		<< "Brute" << space
		<< "PKd" << space
		<< "AKd" << space
		<< "ABd" << space
		<< logNewLine;

	//test(8, 10000, 4, fixed, maxRelativeError);

	std::vector<integer> dimensionSet;
	dimensionSet.push_back(2);
	dimensionSet.push_back(4);
	dimensionSet.push_back(8);
	dimensionSet.push_back(16);
	dimensionSet.push_back(32);
	const integer defaultDimension = 8;

	std::vector<integer> pointCountSet;
	pointCountSet.push_back(2500);
	pointCountSet.push_back(5000);
	pointCountSet.push_back(10000);
	pointCountSet.push_back(20000);
	pointCountSet.push_back(40000);
	const integer defaultPointCount = 10000;

	std::vector<integer> neighborCountSet;
	neighborCountSet.push_back(1);
	neighborCountSet.push_back(2);
	neighborCountSet.push_back(4);
	neighborCountSet.push_back(8);
	neighborCountSet.push_back(16);
	const integer defaultNeighborCount = 4;

	// Vary each dimension independently,
	// while fixing all other dimensions.

	for (integer i = 0;i < dimensionSet.size();++i)
	{
		test(
			dimensionSet[i],
			defaultPointCount,
			defaultNeighborCount,
			fixed, maxRelativeError);
	}

	for (integer i = 0;i < pointCountSet.size();++i)
	{
		test(
			defaultDimension,
			pointCountSet[i],
			defaultNeighborCount,
			fixed, maxRelativeError);
	}
	for (integer i = 0;i < neighborCountSet.size();++i)
	{
		test(
			defaultDimension,
			defaultPointCount,
			neighborCountSet[i], 
			fixed, maxRelativeError);
	}
}

void estimation()
{
	deviceSystem().initialize();

	timLog().addObserver(LogObserverPtr(new StreamLogObserver(&std::cout)));
	timLog().addObserver(LogObserverPtr(new FileLogObserver("timlog.txt")));

	timings();

	deviceSystem().deInitialize();
}

