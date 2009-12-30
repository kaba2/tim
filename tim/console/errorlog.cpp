#include "tim/console/errorlog.h"

namespace Tim
{

	void ErrorLog::report(integer line, const std::string& text)
	{
		//errorMap_.insert(std::make_pair(line, text));
		std::cerr << "Error at line " << line << ": " << text << std::endl;
		exit(1);
	}

	const ErrorLog::Container& ErrorLog::map() const
	{
		return errorMap_;
	}

	std::ostream& operator<<(
		std::ostream& stream, 
		const ErrorLog& errorLog)
	{
		if (!errorLog.map().empty())
		{
			stream << "Errors occurred." << std::endl;

			ErrorLog::ConstIterator iter = errorLog.map().begin();
			const ErrorLog::ConstIterator iterEnd = errorLog.map().end();
			while(iter != iterEnd)
			{
				stream << "Line " << iter->first << ": " << iter->second << std::endl;
				++iter;
			}
		}

		return stream;
	}

	void reportError(const std::string& text)
	{
		std::cerr << text << std::endl;
	}

}
