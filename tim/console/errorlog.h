// Description: ErrorLog class
// Detail: Stores the errors emitted by the interpreter.
// Documentation: tim_console_cpp.txt

#ifndef TIM_ERRORLOG_H
#define TIM_ERRORLOG_H

#include <tim/core/mytypes.h>

#include <string>
#include <map>
#include <iostream>

namespace Tim
{

	class ErrorLog
	{
	public:
		typedef std::multimap<integer, std::string> Container;
		typedef Container::const_iterator ConstIterator;

		ErrorLog();

		void report(const std::string& text);
		void report(integer line, const std::string& text);

		const Container& map() const;
		
		void setLine(integer line);
		integer line() const;

	private:
		 Container errorMap_;
		 integer line_;
	};

	ErrorLog& errorLog();

	std::ostream& operator<<(std::ostream& stream, const ErrorLog& errorLog);

	void reportError(const std::string& text);


}

#endif
