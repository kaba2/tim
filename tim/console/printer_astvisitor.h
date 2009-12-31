// Description: Printer_AstVisitor class
// Detail: Visualizes the program's AST.
// Documentation: tim_console_cpp.txt

#ifndef TIM_PRINTER_ASTVISITOR_H
#define TIM_PRINTER_ASTVISITOR_H

#include "tim/console/astvisitor.h"

#include <iostream>

namespace Tim
{

	class Printer_AstVisitor
		: public AstVisitor
	{
	public:
		explicit Printer_AstVisitor(std::ostream& stream);

		virtual void visit(const Program_AstNode& node);
		virtual void visit(const Statement_AstNode& node);
		virtual void visit(const Expression_AstNode& node);
		virtual void visit(const Declaration_AstNode& node);
		virtual void visit(const Print_AstNode& node);
		virtual void visit(const Identifier_AstNode& node);
		virtual void visit(const Integer_AstNode& node);
		virtual void visit(const Real_AstNode& node);
		virtual void visit(const String_AstNode& node);
		virtual void visit(const CellArray_AstNode& node);
		virtual void visit(const RealArray_AstNode& node);
		virtual void visit(const FunctionCall_AstNode& node);

	private:
		template <typename Type>
		void doVisit(const Type& node);

		std::ostream& output();
		void increaseIndentation();
		void decreaseIndentation();

		std::ostream* stream_;
		integer indentation_;
	};


}

#endif
