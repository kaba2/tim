#include "tim/console/interpreter.h"
#include "tim/console/functions.h"
#include "tim/console/errorlog.h"

#include "tim/core/signal_tools.h"

namespace Tim
{

	void Interpreter_AstVisitor::visit(const Declaration_AstNode& node)
	{
		symbolMap_[node.name()] = evaluate(node.expression());
	}

	void Interpreter_AstVisitor::visit(const Print_AstNode& node)
	{
		boost::any value = evaluate(node.expression());

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
		if (const Identifier_AstNode* node = dynamic_cast<const Identifier_AstNode*>(expression))
		{
			SymbolIterator iter = symbolMap_.find(node->name());
			if (iter == symbolMap_.end())
			{
				reportError("Undefined identifier '" + node->name() + "'.");
				return boost::any();
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

			Cell* cellArray = new Cell(width, height);
			for (integer i = 0;i < cellArray->size();++i)
			{
				const std::string& name = identifierArray(i);
				SymbolIterator iter = symbolMap_.find(name);
				if (iter == symbolMap_.end())
				{
					reportError("Undefined identifier '" + name + "'.");
				}
				else
				{
					try
					{
						(*cellArray)(i) = boost::any_cast<SignalPtr>(iter->second);
					}
					catch(const boost::bad_any_cast&)
					{
						reportError("'" + name + "' is not a signal.");
					}
				}
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
			
			return functionCall(node->name(), inputSet);
		}

		const bool thisPlaceIsReached = true;
		ENSURE(!thisPlaceIsReached);
		
		return boost::any();
	}

}
