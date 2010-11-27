#include "estimation.h"

#include <pastel/device/devicesystem.h>

#include <pastel/sys/log_all.h>
#include <pastel/sys/testreport.h>

#include <iostream>

using namespace Pastel;
using namespace Tim;

void estimation();

int estimationMain()
{
	deviceSystem().initialize();

	log().addLogger(LoggerPtr(new Stream_Logger(&std::cout)));
	log().addLogger(LoggerPtr(new File_Logger("log.txt")));

	timTestList().console();
	generateTestReport(timTestReport(), log());

	deviceSystem().deInitialize();

	return 0;
}
