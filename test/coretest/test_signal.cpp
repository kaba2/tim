#include "estimation.h"

#include "tim/core/signal_tools.h"
#include "tim/core/embed.h"

#include <pastel/sys/view_all.h>
#include <pastel/sys/random.h>
#include <pastel/sys/string_tools.h>

using namespace Tim;

namespace
{

	class SignalTest
		: public TestSuite
	{
	public:
		SignalTest()
			: TestSuite(&timTestReport())
		{
		}

		virtual void run()
		{
			testMerge();
			testEmbed();
			testSlice();
		}

		void testSlice()
		{
			SignalPtr xy = 
				SignalPtr(new Signal(5, 3));

			xy->data() |=
				0, 5, 10,
				1, 6, 11,
				2, 7, 12,
				3, 8, 13,
				4, 9, 14;

			SignalPtr x = split(xy, 0, 1);
			SignalPtr xCorrect =
				SignalPtr(new Signal(5, 1));
			xCorrect->data() |=
				0, 1, 2, 3, 4;

			TEST_ENSURE(x->data() == xCorrect->data());

			SignalPtr y = split(xy, 1, 3);
			SignalPtr yCorrect =
				SignalPtr(new Signal(5, 2));
			yCorrect->data() |=
				5, 10,
				6, 11,
				7, 12,
				8, 13,
				9, 14;

			TEST_ENSURE(y->data() == yCorrect->data());
		}

		void testMerge()
		{
			SignalPtr x =
				SignalPtr(new Signal(5, 1));
			x->data() |= 
				1, 2, 3, 4, 5;

			SignalPtr y =
				SignalPtr(new Signal(5, 1));
			y->data() |= 
				6, 7, 8, 9, 10;

			SignalPtr z = merge(x, y);
			SignalPtr zCorrect =
				SignalPtr(new Signal(5, 2));
			zCorrect->data() |= 
				1, 6,
				2, 7, 
				3, 8,
				4, 9,
				5, 10;

			{
				SignalPtr w = merge(z, z);
				SignalPtr wCorrect =
					SignalPtr(new Signal(5, 4));
				wCorrect->data() |= 
					1, 6, 1, 6,
					2, 7, 2, 7,
					3, 8, 3, 8,
					4, 9, 4, 9, 
					5, 10, 5, 10;

				TEST_ENSURE(w->data() == wCorrect->data());

				SignalPtr u = split(wCorrect, 0, 2);

				TEST_ENSURE(u->data() == zCorrect->data());

				SignalPtr v = split(wCorrect, 2, 4);

				TEST_ENSURE(v->data() == zCorrect->data());
			}
			{
				SignalPtr w = merge(z, z, 1);
				SignalPtr wCorrect =
					SignalPtr(new Signal(4, 4));
				wCorrect->data() |= 
					2, 7, 1, 6,
					3, 8, 2, 7,
					4, 9, 3, 8,
					5, 10, 4, 9;

				TEST_ENSURE(w->data() == wCorrect->data());
			}
			{
				SignalPtr w = merge(z, z, 2);
				SignalPtr wCorrect =
					SignalPtr(new Signal(3, 4));
				wCorrect->data() |= 
					3, 8, 1, 6,
					4, 9, 2, 7,
					5, 10, 3, 8;

				TEST_ENSURE(w->data() == wCorrect->data());
			}
		}

		void testEmbed()
		{
			SignalPtr x = 
				SignalPtr(new Signal(15, 1));
			x->data() |= 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14;

			{
				SignalPtr y = delayEmbed(x, 3);
				SignalPtr yCorrect =
					SignalPtr(new Signal(13, 3));
				yCorrect->data() |=
					0, 1, 2,
					1, 2, 3,
					2, 3, 4,
					3, 4, 5,
					4, 5, 6,
					5, 6, 7,
					6, 7, 8,
					7, 8, 9,
					8, 9, 10,
					9, 10, 11,
					10, 11, 12,
					11, 12, 13,
					12, 13, 14;

				TEST_ENSURE(y->data() == yCorrect->data());
			}

			{
				SignalPtr y = delayEmbed(x, 3, 1);
				SignalPtr yCorrect =
					SignalPtr(new Signal(12, 3));
				yCorrect->data() |=
					1, 2, 3,
					2, 3, 4,
					3, 4, 5,
					4, 5, 6,
					5, 6, 7,
					6, 7, 8,
					7, 8, 9, 
					8, 9, 10,
					9, 10, 11,
					10, 11, 12,
					11, 12, 13,
					12, 13, 14;

				TEST_ENSURE(y->data() == yCorrect->data());
			}

			{
				SignalPtr y = delayEmbed(x, 3, 2);
				SignalPtr yCorrect =
					SignalPtr(new Signal(11, 3));
				yCorrect->data() |=
					2, 3, 4,
					3, 4, 5,
					4, 5, 6,
					5, 6, 7,
					6, 7, 8,
					7, 8, 9, 
					8, 9, 10,
					9, 10, 11,
					10, 11, 12,
					11, 12, 13,
					12, 13, 14;
		
				TEST_ENSURE(y->data() == yCorrect->data());
			}

			{
				SignalPtr y = delayEmbed(x, 3, 0, 2);
				SignalPtr yCorrect =
					SignalPtr(new Signal(11, 3));
				yCorrect->data() |=
					0, 2, 4,
					1, 3, 5,
					2, 4, 6,
					3, 5, 7,
					4, 6, 8,
					5, 7, 9,
					6, 8, 10,
					7, 9, 11,
					8, 10, 12,
					9, 11, 13,
					10, 12, 14;
		
				TEST_ENSURE(y->data() == yCorrect->data());
			}

			{
				SignalPtr y = delayEmbed(x, 3, 1, 2);
				SignalPtr yCorrect =
					SignalPtr(new Signal(10, 3));
				yCorrect->data() |=
					1, 3, 5,
					2, 4, 6,
					3, 5, 7,
					4, 6, 8,
					5, 7, 9,
					6, 8, 10,
					7, 9, 11,
					8, 10, 12,
					9, 11, 13,
					10, 12, 14;
		
				TEST_ENSURE(y->data() == yCorrect->data());
			}
		}
	};

	void testSignal()
	{
		SignalTest test;
		test.run();
	}

	void addTest()
	{
		timTestList().add("Signal", testSignal);
	}

	CallFunction run(addTest);

}

