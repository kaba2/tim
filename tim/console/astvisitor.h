// Description: AstVisitor class
// Detail: A visitor base-class for visiting the AST.
// Documentation: tim_console_cpp.txt

#ifndef TIM_ASTVISITOR_H
#define TIM_ASTVISITOR_H

namespace Tim
{

	class Program_AstNode;
	class Statement_AstNode;
	class Expression_AstNode;
	class Declaration_AstNode;
	class Print_AstNode;
	class Identifier_AstNode;
	class Integer_AstNode;
	class Real_AstNode;
	class String_AstNode;
	class CellArray_AstNode;
	class RealArray_AstNode;
	class FunctionCall_AstNode;

	class AstVisitor
	{
	public:
		virtual ~AstVisitor()
		{
		}

		virtual void visit(const Program_AstNode& node)
		{
		}

		virtual void visit(const Statement_AstNode& node)
		{
		}

		virtual void visit(const Expression_AstNode& node)
		{
		}

		virtual void visit(const Declaration_AstNode& node)
		{
		}

		virtual void visit(const Print_AstNode& node)
		{
		}

		virtual void visit(const Identifier_AstNode& node)
		{
		}

		virtual void visit(const Integer_AstNode& node)
		{
		}

		virtual void visit(const Real_AstNode& node)
		{
		}

		virtual void visit(const String_AstNode& node)
		{
		}

		virtual void visit(const CellArray_AstNode& node)
		{
		}

		virtual void visit(const RealArray_AstNode& node)
		{
		}

		virtual void visit(const FunctionCall_AstNode& node)
		{
		}
	};

}

#endif
