#include "tim/console/interpreter_astvisitor.h"
#include "tim/console/functions.h"
#include "tim/console/errorlog.h"

#include "tim/core/signal_tools.h"

#include <pastel/sys/string_algorithms.h>

namespace Tim
{

	Interpreter_AstVisitor::Interpreter_AstVisitor()
		: symbolMap_()
	{
		symbolMap_.insert(std::make_pair("gaussian", 
			boost::any(generateGaussian(5000, 10))));
	}

	void Interpreter_AstVisitor::visit(const Program_AstNode& node)
	{
		const StatementSet& statementSet = node.statementSet();
		const integer statements = statementSet.size();
		for (integer i = 0;i < statements;++i)
		{
			visit(*statementSet[i]);
		}
	}

	void Interpreter_AstVisitor::visit(const Statement_AstNode& node)
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

	void Interpreter_AstVisitor::visit(const Declaration_AstNode& node)
	{
		symbolMap_[node.name()] = evaluate(node.expression());
	}

	void Interpreter_AstVisitor::visit(const Print_AstNode& node)
	{
		boost::any value = evaluate(node.expression());

		try
		{
			CellPtr cellArray = boost::any_cast<CellPtr>(value);
			std::cout << "{";

			const integer trials = cellArray->width();
			const integer series = cellArray->height();
			for (integer y = 0;y < series;++y)
			{
				if (y > 0)
				{
					std::cout << ";" << std::endl;
				}
				for (integer x = 0;x < trials;++x)
				{
					if (x > 0)
					{
						std::cout << ", ";
					}

					const SignalPtr signal = (*cellArray)(x, y);
					std::cout << "[" << signal->dimension() << " x " 
						<< signal->samples() + signal->t() << " signal]";
					//std::cout << *signal << std::endl;
				}
			}
			std::cout << "}" << std::endl;
		}
		catch(const boost::bad_any_cast&)
		{
		}

		try
		{
			SignalPtr signal = boost::any_cast<SignalPtr>(value);
			std::cout << *signal << std::endl;
		}
		catch(const boost::bad_any_cast&)
		{
		}

		try
		{
			integer k = boost::any_cast<integer>(value);
			std::cout << k << std::endl;
		}
		catch(const boost::bad_any_cast&)
		{
		}

		try
		{
			real k = boost::any_cast<real>(value);
			std::cout << k << std::endl;
		}
		catch(const boost::bad_any_cast&)
		{
		}

		try
		{
			std::string text = boost::any_cast<std::string>(value);
			std::cout << text << std::endl;
		}
		catch(const boost::bad_any_cast&)
		{
		}
	}

	// Private

	boost::any Interpreter_AstVisitor::evaluate(const Expression_AstNode* expression)
	{
		errorLog().setLine(expression->line());

		if (const Identifier_AstNode* node = dynamic_cast<const Identifier_AstNode*>(expression))
		{
			SymbolIterator iter = symbolMap_.find(node->name());
			if (iter == symbolMap_.end())
			{
				reportError("Undefined identifier '" + node->name() + "'.");
				throw Interpreter_Exception();
			}

			return iter->second;
		}
		if (const Real_AstNode* node = dynamic_cast<const Real_AstNode*>(expression))
		{
			return boost::any(node->number());
		}
		if (const Integer_AstNode* node = dynamic_cast<const Integer_AstNode*>(expression))
		{
			return boost::any(node->number());
		}
		if (const String_AstNode* node = dynamic_cast<const String_AstNode*>(expression))
		{
			return boost::any(node->text());
		}
		if (const RealArray_AstNode* node = dynamic_cast<const RealArray_AstNode*>(expression))
		{
			return boost::any(node->signal());
		}
		if (const CellArray_AstNode* node = dynamic_cast<const CellArray_AstNode*>(expression))
		{
			const Array<std::string>& identifierArray = node->cellArray();
			const integer width = identifierArray.width();
			const integer height = identifierArray.height();

			bool errorsFound = false;			
			CellPtr cellArray(new Cell(width, height));
			for (integer y = 0;y < height;++y)
			{
				for (integer x = 0;x < width;++x)
				{
					const std::string& name = identifierArray(x, y);
					SymbolIterator iter = symbolMap_.find(name);
					if (iter == symbolMap_.end())
					{
						reportError("Undefined identifier '" + name + "'" +
							" at element (" + integerToString(y) + ", " + integerToString(x) + ").");
						errorsFound = true;
					}
					else
					{
						try
						{
							(*cellArray)(x, y) = boost::any_cast<SignalPtr>(iter->second);
						}
						catch(const boost::bad_any_cast&)
						{
							reportError("'" + name + "'" + 
								" at element (" + integerToString(y) + ", " + integerToString(x) + 
								") is not a signal.");
							errorsFound = true;
						}
					}
				}
			}
			if (errorsFound)
			{
				throw Interpreter_Exception();
			}
			return boost::any(cellArray);
		}
		if (const FunctionCall_AstNode* node = dynamic_cast<const FunctionCall_AstNode*>(expression))
		{
			AnySet inputSet;

			const std::vector<Expression_AstNode*>& expressionSet = node->expressionSet();
			for (integer i = 0;i < expressionSet.size();++i)
			{
				inputSet.push_back(evaluate(expressionSet[i]));
			}
			
			try
			{
				return functionCall(node->name(), inputSet);
			}
			catch(const FunctionCall_Exception&)
			{
				throw Interpreter_Exception();
			}
		}

		const bool thisPlaceIsReached = true;
		ENSURE(!thisPlaceIsReached);
		
		return boost::any();
	}

}
