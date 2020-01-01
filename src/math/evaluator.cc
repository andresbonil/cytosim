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
class Evaluator
{
public:
    
    /// variable names must be a single letter
    typedef std::pair<char, real> variable;

    /// list of variables
    typedef std::initializer_list<variable> variable_list;
    
    static void skip_space(char const*& str)
    {
        while ( isspace(*str) )
            ++str;
    }

private:
    
    /// list of variables
    variable_list variables_;
    
    real value_(char const*& str)
    {
        for ( variable const& v : variables_ )
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
    
    real number_(char const*& str)
    {
        errno = 0;
        char * end = nullptr;
        real d = strtod(str, &end);
        if ( errno )
            throw InvalidSyntax("Unexpected syntax");
        str = end;
        return d;
    }
    
    real factor_(char const*& str)
    {
        skip_space(str);
        if ( '0' <= *str && *str <= '9' )
            return number_(str);
        else if ( *str == '(' )
        {
            ++str; // '('
            real result = expression_(str);
            if ( *str != ')' )
                throw InvalidSyntax("Unexpected syntax");
            ++str; // ')'
            return result;
        }
        else if ( *str == '-' )
        {
            ++str;
            return -factor_(str);
        }
        else
            return value_(str);
    }
    
    real term_(char const*& str)
    {
        real result = factor_(str);
        while ( 1 )
        {
            skip_space(str);
            char c = *str;
            if ( c == '*' )
                result *= factor_(++str);
            else if ( c == '/' )
                result /= factor_(++str);
            else if ( c == '^' )
                result = pow(result, factor_(++str));
            else
                return result;
        }
    }
    
    real expression_(char const*& str)
    {
        real result = term_(str);
        while ( 1 )
        {
            skip_space(str);
            char c = *str;
            if ( c == '+' )
                result += term_(++str);
            else if ( c == '-' )
                result -= term_(++str);
            else
                return result;
        }
    }
    
public:
    
    Evaluator(variable_list const& v) : variables_(v)
    {
    }
    
    real value(char const*& str)
    {
        //std::clog << "evaluate (" << str << ")\n";
        return expression_(str);
    }
    
    bool inequality(char const*& str)
    {
        //std::clog << "inequality (" << str << ")\n";
        real a = expression_(str);
        skip_space(str);
        char op = *str++;
        if ( op != '<' && op != '>' )
            throw InvalidSyntax("Unexpected syntax");
        if ( *str == '=' )
        {
            ++str;
            real b = expression_(str);
            if ( op == '<' ) return ( a <= b );
            if ( op == '>' ) return ( a >= b );
        }
        real b = expression_(str);
        if ( op == '<' ) return ( a < b );
        if ( op == '>' ) return ( a > b );
        return false;
    }
};
