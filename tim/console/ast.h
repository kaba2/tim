// Description: Abstract Syntax Tree for the program
// Documentation: tim_console_cpp.txt

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

	class Statement_AstNode
		: public AstNodeBase<Statement_AstNode, AstNode>
	{
	};

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

		~Program_AstNode()
		{
			const integer statements = statementSet_->size();
			for (integer i = 0;i < statements;++i)
			{
				delete (*statementSet_)[i];
			}
			delete statementSet_;
		}

		const StatementSet& statementSet() const
		{
			return *statementSet_;
		}

	private:
		const StatementSet* statementSet_;
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

		~Declaration_AstNode()
		{
			delete expression_;
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

		~Print_AstNode()
		{
			delete expression_;
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

		~CellArray_AstNode()
		{
			delete data_;
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

		~FunctionCall_AstNode()
		{
			const integer expressions = expressionSet_->size();
			for (integer i = 0;i < expressions;++i)
			{
				delete (*expressionSet_)[i];
			}
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

}

#include "tim/console/astvisitor.h"

#endif
