#include "estimation.h"

#include "tim/core/signal_tools.h"
#include "tim/core/delay_embed.h"

#include <pastel/sys/view.h>
#include <pastel/sys/random.h>
#include <pastel/sys/string/string_algorithms.h>

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
			Signal xy = 
				Signal(new Signal(5, 3));

			xy->data() |=
				0, 5, 10,
				1, 6, 11,
				2, 7, 12,
				3, 8, 13,
				4, 9, 14;

			Signal x = split(xy, 0, 1);
			Signal xCorrect =
				Signal(new Signal(5, 1));
			xCorrect->data() |=
				0, 1, 2, 3, 4;

			TEST_ENSURE(x->data() == xCorrect->data());

			Signal y = split(xy, 1, 3);
			Signal yCorrect =
				Signal(new Signal(5, 2));
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
			Signal x =
				Signal(new Signal(5, 1));
			x->data() |= 
				1, 2, 3, 4, 5;

			Signal y =
				Signal(new Signal(5, 1));
			y->data() |= 
				6, 7, 8, 9, 10;

			Signal z = merge(x, y);
			Signal zCorrect =
				Signal(new Signal(5, 2));
			zCorrect->data() |= 
				1, 6,
				2, 7, 
				3, 8,
				4, 9,
				5, 10;
			TEST_ENSURE(z->data() == zCorrect->data());
			TEST_ENSURE_OP(z->t(), ==, 0);

			{
				Signal w = merge(z, z);
				Signal wCorrect =
					Signal(new Signal(5, 4));
				wCorrect->data() |= 
					1, 6, 1, 6,
					2, 7, 2, 7,
					3, 8, 3, 8,
					4, 9, 4, 9, 
					5, 10, 5, 10;

				TEST_ENSURE(w->data() == wCorrect->data());
				TEST_ENSURE_OP(w->t(), ==, 0);

				Signal u = split(wCorrect, 0, 2);

				TEST_ENSURE(u->data() == zCorrect->data());
				TEST_ENSURE_OP(u->t(), ==, 0);

				Signal v = split(wCorrect, 2, 4);

				TEST_ENSURE(v->data() == zCorrect->data());
				TEST_ENSURE_OP(v->t(), ==, 0);
			}
			{
				Signal w = merge(z, z, 0, 1);
				Signal wCorrect =
					Signal(new Signal(4, 4));
				wCorrect->data() |= 
					2, 7, 1, 6,
					3, 8, 2, 7,
					4, 9, 3, 8,
					5, 10, 4, 9;

				TEST_ENSURE(w->data() == wCorrect->data());
				TEST_ENSURE_OP(w->t(), ==, 1);
			}
			{
				Signal w = merge(z, z, 0, 2);
				Signal wCorrect =
					Signal(new Signal(3, 4));
				wCorrect->data() |= 
					3, 8, 1, 6,
					4, 9, 2, 7,
					5, 10, 3, 8;

				TEST_ENSURE(w->data() == wCorrect->data());
				TEST_ENSURE_OP(w->t(), ==, 2);
			}
			{
				Signal z2 = Signal(new Signal(*z));
				z2->setT(2);

				Signal w = merge(z, z2);
				Signal wCorrect =
					Signal(new Signal(3, 4));
				wCorrect->data() |= 
					3, 8, 1, 6,
					4, 9, 2, 7,
					5, 10, 3, 8;

				TEST_ENSURE(w->data() == wCorrect->data());
				TEST_ENSURE_OP(w->t(), ==, 2);
			}
		}

		void testEmbed()
		{
			Signal x = 
				Signal(new Signal(15, 1));
			x->data() |= 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14;

			{
				Signal y = delayEmbed(x, 3);
				Signal yCorrect =
					Signal(new Signal(13, 3));
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
				TEST_ENSURE_OP(y->t(), ==, 2);
			}

			{
				Signal y = delayEmbed(x, 3, 2);
				Signal yCorrect =
					Signal(new Signal(11, 3));
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
				TEST_ENSURE_OP(y->t(), ==, 4);
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

