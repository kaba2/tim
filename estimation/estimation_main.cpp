#include <pastel/sys/log.h>
#include <pastel/sys/streamlogobserver.h>
#include <pastel/sys/filelogobserver.h>

#include <ANN/ANN.h>

#include <iostream>

using namespace Pastel;

void estimation();

int main()
{
	//log().addObserver(LogObserverPtr(new StreamLogObserver(&std::cout)));
	log().addObserver(LogObserverPtr(new FileLogObserver("log.txt")));

	estimation();

	return 0;
}
