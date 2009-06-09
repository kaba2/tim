#include "estimation.h"

#include <pastel/device/devicesystem.h>

#include <pastel/sys/log_all.h>

#include <iostream>

using namespace Pastel;
using namespace Tim;

void estimation();

int estimationMain()
{
	deviceSystem().initialize();

	log().addObserver(LogObserverPtr(new StreamLogObserver(&std::cout)));
	log().addObserver(LogObserverPtr(new FileLogObserver("log.txt")));

	timTestList().console();
	//timTestList().run("nearest_performance");

	deviceSystem().deInitialize();

	return 0;
}
