#include "tim/console/printer_astvisitor.h"

namespace Tim
{

	Printer_AstVisitor::Printer_AstVisitor(
		std::ostream& stream)
		: stream_(&stream)
		, indentation_(0)
	{
	}

	void Printer_AstVisitor::visit(const Program_AstNode& node)
	{
		const StatementSet& statementSet = node.statementSet();
		const integer statements = statementSet.size();
		for (integer i = 0;i < statements;++i)
		{
			visit(*statementSet[i]);
		}
	}

	void Printer_AstVisitor::visit(const Statement_AstNode& node)
	{
		if (const Declaration_AstNode* declaration = 
			dynamic_cast<const Declaration_AstNode*>(&node))
		{
			visit(*declaration);
		}
		if (const Print_AstNode* printExpression = 
			dynamic_cast<const Print_AstNode*>(&node))
		{
			visit(*printExpression);
		}
	}

	void Printer_AstVisitor::visit(const Expression_AstNode& node)
	{
		if (const Identifier_AstNode* identifier = 
			dynamic_cast<const Identifier_AstNode*>(&node))
		{
			visit(*identifier);
		}
		if (const Integer_AstNode* integerValue = 
			dynamic_cast<const Integer_AstNode*>(&node))
		{
			visit(*integerValue);
		}
		if (const Real_AstNode* realValue = 
			dynamic_cast<const Real_AstNode*>(&node))
		{
			visit(*realValue);
		}
		if (const String_AstNode* text = 
			dynamic_cast<const String_AstNode*>(&node))
		{
			visit(*text);
		}
		if (const CellArray_AstNode* cellArray = 
			dynamic_cast<const CellArray_AstNode*>(&node))
		{
			visit(*cellArray);
		}
		if (const RealArray_AstNode* realArray = 
			dynamic_cast<const RealArray_AstNode*>(&node))
		{
			visit(*realArray);
		}
		if (const FunctionCall_AstNode* functionCall = 
			dynamic_cast<const FunctionCall_AstNode*>(&node))
		{
			visit(*functionCall);
		}
	}

	void Printer_AstVisitor::visit(const Declaration_AstNode& node)
	{
		output() << "Declaration: " << node.name() << std::endl;
		doVisit(*node.expression());
	}

	void Printer_AstVisitor::visit(const Print_AstNode& node)
	{
		output() << "Print" << std::endl;
		doVisit(*node.expression());

	}

	void Printer_AstVisitor::visit(const Identifier_AstNode& node)
	{
		output() << "Identifier: " << node.name() << std::endl;
	}

	void Printer_AstVisitor::visit(const Integer_AstNode& node)
	{
		output() << "Integer: " << node.number() << std::endl;
	}

	void Printer_AstVisitor::visit(const Real_AstNode& node)
	{
		output() << "Real: " << node.number() << std::endl;
	}

	void Printer_AstVisitor::visit(const String_AstNode& node)
	{
		output() << "String: " << node.text() << std::endl;
	}

	void Printer_AstVisitor::visit(const CellArray_AstNode& node)
	{
		output() << "CellArray" << std::endl;
	}

	void Printer_AstVisitor::visit(const RealArray_AstNode& node)
	{
		output() << "RealArray" << std::endl;
	}

	void Printer_AstVisitor::visit(const FunctionCall_AstNode& node)
	{
		output() << "FunctionCall: " << node.name() << std::endl;
	}

	// Private

	template <typename Type>
	void Printer_AstVisitor::doVisit(const Type& node)
	{
		increaseIndentation();
		visit(node);
		decreaseIndentation();
	}

	std::ostream& Printer_AstVisitor::output()
	{
		for (integer i = 0;i < indentation_;++i)
		{
			(*stream_) << "  ";
		}
		return *stream_;
	}

	void Printer_AstVisitor::increaseIndentation()
	{
		++indentation_;
	}

	void Printer_AstVisitor::decreaseIndentation()
	{
		--indentation_;
	}

}
