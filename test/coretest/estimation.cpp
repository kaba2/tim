#include "estimation.h"

#include <pastel/sys/logging.h>
#include <pastel/sys/testreport.h>

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

	setInvariantFailureAction(
		InvariantFailureAction::Throw);

	timTestList().console();
	generateTestReport(timTestReport(), log());

	return 0;
}
