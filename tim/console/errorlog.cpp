#include "tim/console/errorlog.h"

namespace Tim
{

	ErrorLog::ErrorLog()
		: line_(0)
	{
	}

	void ErrorLog::report(const std::string& text)
	{
		errorMap_.insert(std::make_pair(line_, text));
	}

	void ErrorLog::report(integer line, const std::string& text)
	{
		errorMap_.insert(std::make_pair(line, text));
	}

	const ErrorLog::Container& ErrorLog::map() const
	{
		return errorMap_;
	}

	void ErrorLog::setLine(integer line)
	{
		line_ = line;
	}

	integer ErrorLog::line() const
	{
		return line_;
	}

	ErrorLog& errorLog()
	{
		static ErrorLog theErrorLog;
		return theErrorLog;
	}

	std::ostream& operator<<(
		std::ostream& stream, 
		const ErrorLog& errorLog)
	{
		ErrorLog::ConstIterator iter = errorLog.map().begin();
		const ErrorLog::ConstIterator iterEnd = errorLog.map().end();
		while(iter != iterEnd)
		{
			stream << "Line " << iter->first << ": ";
			stream << iter->second << std::endl;
			++iter;
		}

		return stream;
	}

	void reportError(const std::string& text)
	{
		errorLog().report(text);
	}

}
