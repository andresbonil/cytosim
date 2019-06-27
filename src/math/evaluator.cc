// Cytosim was created by Francois Nedelec. Copyright 2019 Cambridge University.
//
//  evaluator.cc
//
//  Created by Francois Nedelec on 08/02/2019.
//  Copyright 2019 Cambridge University. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

/// a minimal math expression evaluator
/**
 This can evaluate boolean expressions like `X^2 + (Y-3)^3 < 4'
 To support more fancy expressions, we could link here `tinyexpr`:
 https://codeplea.com/tinyexpr
*/
namespace Evaluator
{
    /// variable names must be a single letter
    typedef std::pair<char, real> variable;
    
    /// list of variables
    typedef std::initializer_list<variable> variable_list;
    
    // defined below
    real expression(char const*& str, variable_list const& vars);

    void skip_space(char const*& str)
    {
        while ( isspace(*str) )
            ++str;
    }
    
    real value(char const*& str, variable_list const& vars)
    {
        for ( variable const& v : vars )
        {
            if ( *str == v.first || toupper(*str) == v.first )
            {
                ++str;
                return v.second;
            }
        }
        throw InvalidSyntax("Unknown variable '"+std::string(1,*str)+"'");
        return 0;
    }
    
    real number(char const*& str)
    {
        errno = 0;
        char * end = nullptr;
        real d = strtod(str, &end);
        if ( errno )
            throw InvalidSyntax("Unexpected syntax");
        str = end;
        return d;
    }
    
    real factor(char const*& str, variable_list const& vars)
    {
        skip_space(str);
        if ( '0' <= *str && *str <= '9' )
            return number(str);
        else if ( *str == '(' )
        {
            ++str; // '('
            real result = expression(str, vars);
            if ( *str != ')' )
                throw InvalidSyntax("Unexpected syntax");
            ++str; // ')'
            return result;
        }
        else if ( *str == '-' )
        {
            ++str;
            return -factor(str, vars);
        }
        else
            return value(str, vars);
    }
    
    real term(char const*& str, variable_list const& vars)
    {
        real result = factor(str, vars);
        while ( 1 )
        {
            skip_space(str);
            char c = *str;
            if ( c == '*' )
                result *= factor(++str, vars);
            else if ( c == '/' )
                result /= factor(++str, vars);
            else if ( c == '^' )
                result = pow(result, factor(++str, vars));
            else
                return result;
        }
    }
    
    real expression(char const*& str, variable_list const& vars)
    {
        real result = term(str, vars);
        while ( 1 )
        {
            skip_space(str);
            char c = *str;
            if ( c == '+' )
                result += term(++str, vars);
            else if ( c == '-' )
                result -= term(++str, vars);
            else
                return result;
        }
    }
    
    real evaluate(char const*& str, variable_list const& vars)
    {
        //std::clog << "evaluate (" << str << ")\n";
        return expression(str, vars);
    }
    
    bool inequality(char const*& str, variable_list const& vars)
    {
        //std::clog << "inequality (" << str << ")\n";
        real a = expression(str, vars);
        skip_space(str);
        char op = *str++;
        if ( op != '<' && op != '>' )
            throw InvalidSyntax("Unexpected syntax");
        if ( *str == '=' )
        {
            ++str;
            real b = expression(str, vars);
            if ( op == '<' ) return ( a <= b );
            if ( op == '>' ) return ( a >= b );
        }
        real b = expression(str, vars);
        if ( op == '<' ) return ( a < b );
        if ( op == '>' ) return ( a > b );
        return false;
    }
}
