// Description: Interpreter_AstVisitor class
// Detail: Interprets the program stored in the AST.
// Documentation: tim_console_cpp.txt

#ifndef TIM_INTERPRETER_ASTVISITOR_H
#define TIM_INTERPRETER_ASTVISITOR_H

#include "tim/console/astvisitor.h"

#include <map>
#include <string>

#include <boost/any.hpp>

namespace Tim
{

	class Interpreter_AstVisitor
		: public AstVisitor
	{
	public:
		Interpreter_AstVisitor();

		virtual void visit(const Program_AstNode& node);
		virtual void visit(const Statement_AstNode& node);
		virtual void visit(const Declaration_AstNode& node);
		virtual void visit(const Print_AstNode& node);

	private:
		typedef std::map<std::string, boost::any> SymbolMap;
		typedef SymbolMap::const_iterator SymbolIterator;

		boost::any evaluate(const Expression_AstNode* expression);

		SymbolMap symbolMap_;
	};

	class Interpreter_Exception
	{
	};

}

#endif
