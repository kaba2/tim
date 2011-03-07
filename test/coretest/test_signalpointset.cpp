#include "estimation.h"

#include "tim/core/signalpointset.h"
#include "tim/core/signal_tools.h"
#include "tim/core/embed.h"

#include <pastel/sys/view_all.h>
#include <pastel/sys/random.h>
#include <pastel/sys/string_algorithms.h>

using namespace Tim;

namespace
{

	class SignalPointSetTest
		: public TestSuite
	{
	public:
		SignalPointSetTest()
			: TestSuite(&timTestReport())
		{
		}

		virtual void run()
		{
			testBasic();
			testBasic2();
		}

		void testBasic()
		{
			SignalPtr xy = 
				SignalPtr(new Signal(5, 3));

			xy->data() |=
				0, 5, 10,
				1, 6, 11,
				2, 7, 12,
				3, 8, 13,
				4, 9, 14;

			SignalPtr signalSet[] = {xy};

			SignalPointSet pointSet(range(signalSet));

			TEST_ENSURE_OP(pointSet.windowBegin() + pointSet.samples(), ==, pointSet.windowEnd());
			TEST_ENSURE_OP(std::distance(pointSet.begin(), pointSet.end()), ==, xy->samples());
			TEST_ENSURE_OP(pointSet.kdTree().points(), ==, xy->samples());
			TEST_ENSURE_OP(pointSet.samples(), ==, xy->samples());
			TEST_ENSURE_OP(pointSet.dimension(), ==, xy->dimension());

			changeTimeWindow(pointSet, -10, 10);
			changeTimeWindow(pointSet, 0, 10);
			changeTimeWindow(pointSet, 1, 10);
			changeTimeWindow(pointSet, 1, 9);
			changeTimeWindow(pointSet, 1, 4);
			changeTimeWindow(pointSet, 1, 3);
			changeTimeWindow(pointSet, 1, 2);
			changeTimeWindow(pointSet, 1, 1);
			changeTimeWindow(pointSet, 2, 2);
			changeTimeWindow(pointSet, 2, 5);
			changeTimeWindow(pointSet, 0, 5);
			changeTimeWindow(pointSet, 0, 1);
			changeTimeWindow(pointSet, 0, 10);
			changeTimeWindow(pointSet, -10, 10);
			changeTimeWindow(pointSet, -20, 10);
			changeTimeWindow(pointSet, -20, 20);
			changeTimeWindow(pointSet, 2, 3);
			changeTimeWindow(pointSet, 1, 4);
			changeTimeWindow(pointSet, 1, 4);

			TEST_ENSURE((*pointSet.begin())->point() == &(xy->data()(3 * 1)));
			TEST_ENSURE((*(pointSet.end() - 1))->point() == &(xy->data()(3 * 3)));

			{
				SignalPointSet pointSet(range(signalSet));

				TEST_ENSURE_OP(std::distance(pointSet.begin(), pointSet.end()), ==, xy->samples());
				TEST_ENSURE_OP(pointSet.kdTree().points(), ==, xy->samples());
				TEST_ENSURE_OP(pointSet.samples(), ==, xy->samples());
				TEST_ENSURE_OP(pointSet.dimension(), ==, xy->dimension());
			}
			{
				SignalPointSet pointSet(range(signalSet), 1, 2);

				TEST_ENSURE_OP(std::distance(pointSet.begin(), pointSet.end()), ==, xy->samples());
				TEST_ENSURE_OP(pointSet.kdTree().points(), ==, xy->samples());
				TEST_ENSURE_OP(pointSet.samples(), ==, xy->samples());
				TEST_ENSURE_OP(pointSet.dimension(), ==, 1);

				TEST_ENSURE((*pointSet.begin())->point() == &(xy->data()(3 * 0 + 1)));
				TEST_ENSURE((*(pointSet.end() - 1))->point() == &(xy->data()(3 * 4 + 1)));
			}
		}

		void testBasic2()
		{
			SignalPtr xy = 
				SignalPtr(new Signal(5, 3, 1));

			xy->data() |=
				0, 5, 10,
				1, 6, 11,
				2, 7, 12,
				3, 8, 13,
				4, 9, 14;

			SignalPtr signalSet[] = {xy};

			SignalPointSet pointSet(range(signalSet));

			TEST_ENSURE_OP(pointSet.windowBegin(), ==, 1);
			TEST_ENSURE_OP(pointSet.windowEnd(), ==, 6);
			TEST_ENSURE_OP(std::distance(pointSet.begin(), pointSet.end()), ==, xy->samples());
			TEST_ENSURE_OP(pointSet.samples(), ==, xy->samples());
			TEST_ENSURE_OP(pointSet.dimension(), ==, xy->dimension());

			changeTimeWindow(pointSet, -10, 10);
			changeTimeWindow(pointSet, 0, 10);
			changeTimeWindow(pointSet, 1, 10);
			changeTimeWindow(pointSet, 1, 9);
			changeTimeWindow(pointSet, 1, 4);
			changeTimeWindow(pointSet, 1, 3);
			changeTimeWindow(pointSet, 1, 2);
			changeTimeWindow(pointSet, 1, 1);
			changeTimeWindow(pointSet, 2, 2);
			changeTimeWindow(pointSet, 2, 5);
			changeTimeWindow(pointSet, 0, 5);
			changeTimeWindow(pointSet, 0, 1);
			changeTimeWindow(pointSet, 0, 10);
			changeTimeWindow(pointSet, -10, 10);
			changeTimeWindow(pointSet, -20, 10);
			changeTimeWindow(pointSet, -20, 20);
			changeTimeWindow(pointSet, 2, 3);
			changeTimeWindow(pointSet, 1, 4);
			changeTimeWindow(pointSet, 1, 4);

			TEST_ENSURE((*pointSet.begin())->point() == &(xy->data()(3 * 0)));
			TEST_ENSURE((*(pointSet.end() - 1))->point() == &(xy->data()(3 * 2)));
		}

		void changeTimeWindow(SignalPointSet& pointSet, integer begin, integer end)
		{
			pointSet.setTimeWindow(begin, end);
			const integer tWidth = pointSet.windowEnd() - pointSet.windowBegin();
			TEST_ENSURE_OP(std::distance(pointSet.begin(), pointSet.end()), ==, tWidth);
			TEST_ENSURE_OP(pointSet.kdTree().points(), ==, tWidth);
		}
	};

	void testSignalPointSet()
	{
		SignalPointSetTest test;
		test.run();
	}

	void addTest()
	{
		timTestList().add("SignalPointSet", testSignalPointSet);
	}

	CallFunction run(addTest);

}

