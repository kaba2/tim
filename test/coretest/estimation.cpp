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
	log().addLogger(LoggerPtr(new Stream_Logger(&std::cout)));
	log().addLogger(LoggerPtr(new File_Logger("log.txt")));

	setInvariantFailureAction(
		InvariantFailureAction::Throw);

	deviceSystem().initialize();

	timTestList().console();
	generateTestReport(timTestReport(), log());

	deviceSystem().deInitialize();

	return 0;
}
