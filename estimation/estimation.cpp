#include <ANN/ANN.h>

#include <boost/random.hpp>

#include <pastel/device/devicesystem.h>
#include <pastel/device/timer.h>
#include <pastel/device/mutex.h>
#include <pastel/device/thread.h>

#include <pastel/sys/vector_tools.h>
#include <pastel/sys/log.h>
#include <pastel/sys/random.h>

#include <pastel/math/uniformsampling.h>

#include <cmath>

using namespace std;
using namespace Pastel;

template <int N, typename Real>
void generateBallSet(integer points,
					 std::vector<Vector<N, Real> >& pointSet)
{
	generateBallSet(point, N, pointSet);
}

template <int N, typename Real>
void generateBallSet(integer points,
					 integer dimension,
					 std::vector<Vector<N, Real> >& pointSet)
{
	std::vector<Vector<N, Real> > result;
	result.reserve(points);

	for (integer i = 0;i < points;++i)
	{
		result.push_back(randomVectorBall<N, Real>(dimension));
	}

	result.swap(pointSet);
}

class NearestComputation
{
public:
	friend class NearestThread;

	NearestComputation(
		const std::vector<Vector<Unbounded, real> >& pointSet,
		integer kNearest,
		real epsilon);

	~NearestComputation();

	void work();

private:
	ANNpointArray		annPointSet_;
	ANNkd_tree*			kdTree_;

	integer points_;
	integer dimension_;
	integer kNearest_;
	real epsilon_;
};

class NearestThread
	: public Thread
{
public:
	explicit NearestThread(
		NearestComputation* computation,
		integer indexBegin,
		integer indexEnd)
		: computation_(computation)
		, indexBegin_(indexBegin)
		, indexEnd_(indexEnd)
		, nearestSet_(0)
		, distanceSet_(0)
	{
		const integer kNearest = computation_->kNearest_;
		nearestSet_ = new ANNidx[kNearest];
		distanceSet_ = new ANNdist[kNearest];
	}

	~NearestThread()
	{
		delete [] nearestSet_;
		delete [] distanceSet_;
	}

	virtual void run()
	{
		ANNkd_tree* kdTree = computation_->kdTree_;
		
		ANNpointArray pointSet = computation_->annPointSet_; 
		const integer kNearest = computation_->kNearest_;
		const integer epsilon = computation_->epsilon_;

		for (integer i = indexBegin_;i < indexEnd_;++i)
		{
			kdTree->annkSearch(	
				pointSet[i],
				kNearest,	
				nearestSet_,				
				distanceSet_,
				epsilon);
		}
	}

private:
	NearestComputation* computation_;
	integer indexBegin_;
	integer indexEnd_;
	ANNidxArray nearestSet_;
	ANNdistArray distanceSet_;
};

typedef CountedPtr<NearestThread> NearestThreadPtr;

NearestComputation::NearestComputation(const std::vector<Vector<Unbounded, real> >& pointSet,
									   integer kNearest,
									   real epsilon)

{
	ENSURE(!pointSet.empty());

	++kNearest;

	points_ = pointSet.size();
	dimension_ = pointSet.front().size();
	epsilon_ = epsilon;
	kNearest_ = kNearest;

	annPointSet_ = new ANNpoint[points_];

	for (integer i = 0;i < points_;++i)
	{
		annPointSet_[i] = (real*)&pointSet[i][0];
	}

	kdTree_ = new ANNkd_tree(
		annPointSet_,		
		points_,			
		dimension_);				
}

void NearestComputation::work()
{
	Timer timer;
	timer.setStart();

	std::vector<NearestThreadPtr> threadSet;

	integer threads = 2;

	// Divide the work evenly to threads.

	integer indexBegin = 0;
	const real indexAdd = points_ / threads;

	for (integer i = 0;i < threads;++i) 
	{
		integer indexEnd = std::min(
			(integer)(indexBegin + indexAdd),
			points_);
		if (i == threads - 1)
		{
			indexEnd = points_;
		}

		threadSet.push_back(
			NearestThreadPtr(new NearestThread(this, indexBegin, indexEnd)));

		log() << "Thread " << i << ": [" << indexBegin 
			<< ", " << indexEnd << "[" << logNewLine;

		indexBegin = indexEnd;
	}

	// Send the threads to work.

	for (integer i = 0;i < threads;++i)
	{
		threadSet[i]->launch();
	}

	// Wait for the threads to finish.

	for (integer i = 0;i < threads;++i)
	{
		threadSet[i]->wait();
	}

	timer.store();

	log() << dimension_ << "-D, " << points_ << " points, " 
		<< kNearest_ << " neighbors." << logNewLine;
	log() << timer.seconds() << " seconds." << logNewLine;
}

NearestComputation::~NearestComputation()
{
	delete [] annPointSet_;
	annPointSet_ = 0;
	delete kdTree_;
	kdTree_ = 0;
}

void test(integer dimension, integer points, integer kNearest, real eps = 0)
{
	std::vector<Vector<Unbounded, real> > pointSet;
	generateBallSet(points, dimension, pointSet);

	NearestComputation nearestComputation(pointSet, kNearest, eps);
	nearestComputation.work();
}

void timings()
{
	const real eps = 0;

	test(5, 5000, 1, eps);
	test(10, 5000, 1, eps);
	test(20, 5000, 1, eps);

	test(5, 10000, 1, eps);
	test(10, 10000, 1, eps);
	test(20, 10000, 1, eps);

	test(5, 20000, 1, eps);
	test(10, 20000, 1, eps);
	test(20, 20000, 1, eps);

	test(5, 10000, 4, eps);
	test(10, 10000, 4, eps);
	test(20, 10000, 4, eps);

	test(5, 10000, 8, eps);
	test(10, 10000, 8, eps);
	test(20, 10000, 8, eps);
}

void estimation()
{
	deviceSystem().initialize();

	timings();

	annClose();					
	deviceSystem().deInitialize();
}

