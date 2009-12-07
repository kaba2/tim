#ifndef TIM_ERROR_H
#define TIM_ERROR_H

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

		void report(integer line, const std::string& text);

		const Container& map() const;	

	private:
		 Container errorMap_;
	};

	std::ostream& operator<<(std::ostream& stream, const ErrorLog& errorLog);


}

#endif
