// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "stream_func.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cctype>


void StreamFunc::clean_stream(std::ostream& os, std::istream& is)
{
    char c = 0;
    while ( is.good() )
    {
        is.get(c);
        
        // terminate the line for new-line and cariage-return
        if ( c == '\r' )
            os << '\n';
        // all type of spaces are substituted
        else if ( isspace(c) )
            os << ' ';
        // non=printable characters are removed
        else if ( isprint(c) )
            os << c;
        else
            std::cerr << "unprintable ascii "<< (int)c << " found\n";
    }
}


void StreamFunc::diff_stream(std::ostream& os, std::istream& val, std::istream& ref)
{
    std::string val_l, ref_l;
    val.seekg(0);
    ref.seekg(0);
    
    while ( val.good() )
    {
        std::getline(val, val_l);
        std::getline(ref, ref_l);
#if ( 0 )
        // print any line containing '{' or '}' 
        bool par = ( std::string::npos != ref_l.find_first_of("{}") );
        if ( val_l != ref_l || par )
#else
        if ( val_l != ref_l )
#endif
        {
            os << val_l << '\n';
            //os << val_l << "  (" << ref_l << ")\n";
        }
    }
}


void StreamFunc::skip_lines(std::ostream& os, std::istream& is, char skip)
{
    std::string line;
    
    while ( is.good() )
    {
        std::getline(is, line);
        if ( line.size() < 1 )
            continue;
        else if ( line[0] != skip )
            os << line << '\n';
    }
}


void StreamFunc::prefix_lines(std::ostream& os, std::istream& is, const char prefix[],
                              char keep, char skip)
{
    std::string line;
    
    while ( is.good() )
    {
        std::getline(is, line);
        if ( line.size() < 1 )
            continue;
        else if ( line[0] == keep )
            os << line << '\n';
        else if ( line[0] == skip )
            ;
        else
            os << prefix << line << '\n';
    }
}


/**
 The alignment of the vertical bar should match the one in PREF
 */
void print_line(std::ostream& os, int num, std::string const& line)
{
    os << std::setw(9) << num << " | " << line << '\n';
}

/**
 The alignment of the vertical bar should match the one in PREF
 */
void print_line(std::ostream& os, const char prefix[], std::string const& line)
{
    os << prefix << " " << line << '\n';
}

/**
 Print the one line extracted from `is` containing `pos` and indicate
 the position `pos` with a arrowhead in a second line below.
 */
void StreamFunc::mark_line(std::ostream& os, std::istream& is, std::streampos pos, const char prefix[])
{
    if ( !is.good() )
        is.clear();
    
    std::streampos isp, sos = is.tellg();
    if ( pos < 0 ) pos = sos;
    is.seekg(0);

    std::string line;
    
    while ( is.good()  &&  is.tellg() <= pos )
    {
        isp = is.tellg();
        std::getline(is, line);
    }
    
    print_line(os, prefix, line);
    is.clear();
    is.seekg(isp);
    
    line.clear();
    int c = 0;
    while ( is.tellg() < pos )
    {
        c = is.get();
        if ( c == EOF )
            break;
        if ( isspace(c) )
            line.push_back((char)c);
        else
            line.push_back(' ');
    }
    line.push_back('^');
    print_line(os, prefix, line);

    is.clear();
    is.seekg(sos);
}


std::string StreamFunc::marked_line(std::istream& is, std::streampos pos, const char prefix[])
{
    std::ostringstream oss;
    mark_line(oss, is, pos, prefix);
    return oss.str();
}


/**
 Output enough lines to cover the area specified by [start, end].
 Each line is printed in full and preceded with a line number
 */
void StreamFunc::print_lines(std::ostream& os, std::istream& is,
                             std::streampos start, std::streampos end)
{
    if ( !is.good() )
        is.clear();
    
    std::streampos isp = is.tellg();
    is.seekg(0);
    std::string line;
    
    unsigned int cnt = 0;
    while ( is.good()  &&  is.tellg() <= start  )
    {
        std::getline(is, line);
        ++cnt;
    }

    print_line(os, cnt, line);
    while ( is.good() &&  is.tellg() < end )
    {
        std::getline(is, line);
        ++cnt;
        if ( !std::all_of(line.begin(),line.end(),isspace) )
            print_line(os, cnt, line);
    }

    is.clear();
    is.seekg(isp);
}


std::string StreamFunc::get_lines(std::istream& is, std::streampos s, std::streampos e)
{
    std::ostringstream oss;
    print_lines(oss, is, s, e);
    return oss.str();
}


unsigned StreamFunc::line_number(std::istream& is, std::streampos pos)
{
    if ( !is.good() )
        is.clear();
    
    std::streampos sos = is.tellg();
    is.seekg(0);
    
    if ( pos == -1 )
        pos = sos;
    
    unsigned cnt = 0;
    std::string line;
    
    while ( is.good()  &&  is.tellg() <= pos )
    {
        std::getline(is, line);
        ++cnt;
    }
    
    is.clear();
    is.seekg(sos);
    return cnt;
}


int StreamFunc::find_and_replace(std::string & src,
                                 std::string const& fnd, std::string const& rep)
{
    int num = 0;
    std::string::size_type fLen = fnd.size();
    std::string::size_type rLen = rep.size();
    std::string::size_type pos = src.find(fnd, 0);
    while ( pos != std::string::npos )
    {
        num++;
        src.replace(pos, fLen, rep);
        pos += rLen;
        pos = src.find(fnd, pos);
    }
    return num;
}


