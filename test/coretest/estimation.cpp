#include "estimation.h"

#include <pastel/sys/logging.h>
#include <pastel/sys/testing/testreport.h>

#include <iostream>

using namespace Pastel;
using namespace Tim;

void estimation();

int estimationMain()
{
	Stream_Logger streamLogger(&std::cout);
	File_Logger fileLogger("log.txt");

	log().addLogger(&streamLogger);
	log().addLogger(&fileLogger);

	timTestList().console();
	generateTestReport(timTestReport(), log());

	return 0;
}
