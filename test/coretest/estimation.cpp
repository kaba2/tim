#include "estimation.h"

#include <pastel/device/devicesystem.h>

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

	deviceSystem().initialize();

	timTestList().console();
	generateTestReport(timTestReport(), log());

	deviceSystem().deInitialize();

	return 0;
}
