#include "tim/console/errorlog.h"

namespace Tim
{

	void ErrorLog::report(integer line, const std::string& text)
	{
		errorMap_.insert(std::make_pair(line, text));
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

}
