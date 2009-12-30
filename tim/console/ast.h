#ifndef TIM_AST_H
#define TIM_AST_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include <pastel/sys/array.h>

#include <vector>

namespace Tim
{

	class AstVisitor;

	class AstNode
	{
	public:
		virtual ~AstNode()
		{
		}

		virtual void accept(AstVisitor& visitor) = 0;
	};

	template <typename Derived, typename Base>
	class AstNodeBase
		: public Base
	{
	public:
		virtual void accept(AstVisitor& visitor)
		{
			visitor.visit((Derived&)*this);
		}
	};

	class Statement_AstNode;
	typedef std::vector<Statement_AstNode*> StatementSet;

	class Program_AstNode
		: public AstNodeBase<Program_AstNode, AstNode>
	{
	public:
		explicit Program_AstNode(
			const StatementSet* statementSet)
			: statementSet_(statementSet)
		{
		}

		const StatementSet& statementSet() const
		{
			return *statementSet_;
		}

	private:
		const StatementSet* statementSet_;
	};

	class Statement_AstNode
		: public AstNodeBase<Statement_AstNode, AstNode>
	{
	};

	// Statements

	class Expression_AstNode
		: public AstNodeBase<Expression_AstNode, Statement_AstNode>
	{
	};

	typedef std::vector<Expression_AstNode*> ExpressionSet;

	class Declaration_AstNode
		: public AstNodeBase<Declaration_AstNode, Statement_AstNode>
	{
	public:
		Declaration_AstNode(
			const std::string& name,
			Expression_AstNode* expression)
			: name_(name)
			, expression_(expression)
		{
		}

		const std::string& name() const
		{
			return name_;
		}

		Expression_AstNode* expression() const
		{
			return expression_;
		}

	private:
		std::string name_;
		Expression_AstNode* expression_;
	};

	class Print_AstNode
		: public AstNodeBase<Print_AstNode, Statement_AstNode>
	{
	public:
		explicit Print_AstNode(
			Expression_AstNode* expression)
			: expression_(expression)
		{
		}

		Expression_AstNode* expression() const
		{
			return expression_;
		}

	private:
		Expression_AstNode* expression_;
	};

	// Expressions

	class Identifier_AstNode
		: public AstNodeBase<Identifier_AstNode, Expression_AstNode>
	{
	public:
		explicit Identifier_AstNode(
			const std::string& name)
			: name_(name)
		{
		}

		const std::string& name() const
		{
			return name_;
		}

	private:
		std::string name_;
	};

	class Integer_AstNode
		: public AstNodeBase<Integer_AstNode, Expression_AstNode>
	{
	public:
		explicit Integer_AstNode(
			integer number)
			: number_(number)
		{
		}

		integer number() const
		{
			return number_;
		}

	private:
		integer number_;
	};

	class Real_AstNode
		: public AstNodeBase<Real_AstNode, Expression_AstNode>
	{
	public:
		explicit Real_AstNode(
			real number)
			: number_(number)
		{
		}

		real number() const
		{
			return number_;
		}

	private:
		real number_;
	};

	class String_AstNode
		: public AstNodeBase<String_AstNode, Expression_AstNode>
	{
	public:
		explicit String_AstNode(
			const std::string& text)
			: text_(text)
		{
		}

		const std::string& text() const
		{
			return text_;
		}

	private:
		std::string text_;
	};

	class CellArray_AstNode
		: public AstNodeBase<CellArray_AstNode, Expression_AstNode>
	{
	public:
		explicit CellArray_AstNode(
			Array<std::string>* data)
			: data_(data)
		{
		}

		const Array<std::string>& cellArray() const
		{
			return *data_;
		}

	private:
		Array<std::string>* data_;
	};

	class RealArray_AstNode
		: public AstNodeBase<RealArray_AstNode, Expression_AstNode>
	{
	public:
		explicit RealArray_AstNode(
			const SignalPtr& signal)
			: signal_(signal)
		{
		}

		const SignalPtr& signal() const
		{
			return signal_;
		}

	private:
		SignalPtr signal_;
	};

	class FunctionCall_AstNode
		: public AstNodeBase<FunctionCall_AstNode, Expression_AstNode>
	{
	public:
		explicit FunctionCall_AstNode(
			const std::string& name,
			const ExpressionSet* expressionSet)
			: name_(name)
			, expressionSet_(expressionSet)
		{
		}

		const std::string& name() const
		{
			return name_;
		}

		const ExpressionSet& expressionSet() const
		{
			return *expressionSet_;
		}

	private:
		std::string name_;
		const ExpressionSet* expressionSet_;
	};

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
