#include "estimation.h"

#include <pastel/sys/log_all.h>

#include <iostream>

using namespace Pastel;
using namespace Tim;


void estimation();

int main()
{
	log().addObserver(LogObserverPtr(new StreamLogObserver(&std::cout)));
	log().addObserver(LogObserverPtr(new FileLogObserver("log.txt")));

	timTestList().console();

	//estimation();

	return 0;
}
