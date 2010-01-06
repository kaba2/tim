#include "tim/console/errorlog.h"

namespace Tim
{

	ErrorLog::ErrorLog()
		: errorMap_()
		, line_(0)
		, nameStack_()
	{
	}

	void ErrorLog::pushNamespace(const std::string& name)
	{
		nameStack_.push_back(name);
	}

	void ErrorLog::popNamespace()
	{
		nameStack_.pop_back();
	}

	std::string ErrorLog::prefix() const
	{
		std::string namePrefix;
		const integer names = nameStack_.size();
		for (integer i = 0;i < names;++i)
		{
			namePrefix += nameStack_[i];
		}
		
		return namePrefix;
	}

	void ErrorLog::report(const std::string& text)
	{
		errorMap_.insert(std::make_pair(line_, prefix() + text));
	}

	void ErrorLog::report(integer line, const std::string& text)
	{
		errorMap_.insert(std::make_pair(line, prefix() + text));
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

	void reportError(integer line, const std::string& text)
	{
		errorLog().report(line, text);
	}

	void reportError(const std::string& text)
	{
		errorLog().report(text);
	}

}
