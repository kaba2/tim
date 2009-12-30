#ifndef TIM_INTERPRETER_H
#define TIM_INTERPRETER_H

#include "tim/console/ast.h"

#include <map>
#include <string>

#include <boost/any.hpp>

namespace Tim
{

	class Interpreter_AstVisitor
		: public AstVisitor
	{
	public:
		virtual void visit(const Program_AstNode& node);
		virtual void visit(const Statement_AstNode& node);
		virtual void visit(const Declaration_AstNode& node);
		virtual void visit(const Print_AstNode& node);

	private:
		boost::any evaluate(const Expression_AstNode* expression);

		typedef std::map<std::string, boost::any> SymbolMap;
		typedef SymbolMap::const_iterator SymbolIterator;
		SymbolMap symbolMap_;
	};

	class Interpreter_Exception
	{
	};

}

#endif
