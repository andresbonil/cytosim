// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef PARSER_H
#define PARSER_H

#include "interface.h"


/// Parser to read and execute Cytosim config files
/**
 This is where the syntax of the config file is defined
 */
class Parser : public Interface
{
private:
    
    /// disabled default constructor
    Parser();
    
    /// control switch to enable command 'set' (creating a property)
    bool      do_set;
    
    /// control switch to enable command 'change' (change a property)
    bool      do_change;
    
    /// control switch to enable command 'new' and 'delete' (object)
    bool      do_new;
    
    /// control switch to enable command 'run' (run simulation)
    bool      do_run;
    
    /// control switch to enable command 'write' (write files)
    bool      do_write;
 
    /// position of stream at the start of current parsing task
    std::streampos spos;
    
    /// print the lines located between `pos` and current position
    void show_lines(std::istream&, std::streampos pos);

public:
    
    /// set the permission of the parser
    Parser(Simul& s, bool allow_set, bool allow_change, bool allow_new, bool allow_run, bool allow_write);
    
    //-------------------------------------------------------------------------------
    
    /// parse command \b set
    void      parse_set(std::istream&);
    
    /// parse command \b change
    void      parse_change(std::istream&);
    
    /// parse command \b new
    void      parse_new(std::istream&);
    
    /// parse command \b delete
    void      parse_delete(std::istream&);
    
    /// parse command \b mark
    void      parse_mark(std::istream&);

    /// parse command \b cut
    void      parse_cut(std::istream&);

    /// parse command \b run
    void      parse_run(std::istream&);
    
    /// parse command \b read
    void      parse_read(std::istream&);
    
    /// parse command \b write
    void      parse_report(std::istream&);
    
    /// parse command \b import
    void      parse_import(std::istream&);
    
    /// parse command \b export
    void      parse_export(std::istream&);
    
    /// parse command \b call
    void      parse_call(std::istream&);
    
    /// parse command \b repeat
    void      parse_repeat(std::istream&);

    /// parse command \b for
    void      parse_for(std::istream&);
    
    /// parse command \b end
    void      parse_end(std::istream&);

    
    /// Parse content of stream
    void      evaluate(std::istream&);

    /// Parse stream, using `msg` to report errors
    void      evaluate(std::istream&, std::string const& msg);

    /// Parse code, using `msg` to report errors
    void      evaluate(std::string const&, std::string const& msg);

    //-------------------------------------------------------------------------------

    /// Open and parse the config file with the given name
    int       readConfig(std::string const& name);

    /// Parse the default config file i.e. calls readConfig(simul.prop->config)
    int       readConfig();

};

#endif

