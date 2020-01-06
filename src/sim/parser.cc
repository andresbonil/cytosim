// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "sim.h"
#include "parser.h"
#include "messages.h"
#include "tokenizer.h"
#include "glossary.h"
#include "filepath.h"
#include "stream_func.h"
#include "simul_prop.h"
#include <fstream>


// Use the second definition to get some verbose reports:
#define VLOG(ARG) ((void) 0)
//#define VLOG(ARG) std::clog << ARG;


//------------------------------------------------------------------------------
/**
 The permission of the parser are:
 - do_set: new Properties can be created
 - do_change: existing Property or Object can be modified
 - do_new: new Object can be created
 - do_run: can perform simulation steps
 - do_write: can write to disc
 .
 */
Parser::Parser(Simul& s, bool ds, bool dc, bool dn, bool dr, bool dw)
: Interface(s), do_set(ds), do_change(dc), do_new(dn), do_run(dr), do_write(dw)
{
    //std::clog << " set " << ds << " change " << dc << " new " << dn << " run " << dr << " write " << dw << "\n";
}


void Parser::show_lines(std::istream& is, std::streampos pos)
{
    std::cerr << "  in\n";
    StreamFunc::print_lines(std::cerr, is, pos, is.tellg());
}


//------------------------------------------------------------------------------
#pragma mark - Parse

/**
 Create a new Property, which is a set of parameters associated with a class.
 
     set CLASS NAME
     {
       PARAMETER = VALUE
       ...
     }
 
 CLASS should be one of the predefined object class (see @ref ObjectGroup).\n
 NAME should be a string, starting with a letter, followed by alphanumeric characters.
 The underscore character is also allowed, as in `fiber_23`\n
 
 It is also possible to use 'set' to change values of an existing Property:
 
     set NAME
     {
       PARAMETER = VALUE
       ...
     }
 
 This is equivalent to:
 
     set NAME { PARAMETER = VALUE }
 
 or:
 
     set NAME PARAMETER { VALUE }
 
 */

void Parser::parse_set(std::istream& is)
{
    std::streampos ipos = is.tellg();
    std::string cat = Tokenizer::get_symbol(is);
    std::string name, para, blok;
    
#ifdef BACKWARD_COMPATIBILITY
    {
        // Read formats anterior to 3.11.2017 ('set hand 2 kinesin')
        unsigned inx = 0;
        Tokenizer::get_integer(is, inx);
    }
#endif
    
    Glossary opt;
    Property * pp = nullptr;

    bool spec = ( is.peek() == ':' );
    
    if ( simul.isPropertyClass(cat) && !spec )
    {
        /* in this form, 'set' defines a new Property
         set CLASS NAME { PARAMETER = VALUE }
         */
        name = Tokenizer::get_symbol(is);
        blok = Tokenizer::get_block(is, '{');
        
        if ( blok.empty() )
            throw InvalidSyntax("syntax error");

        if ( do_set )
        {
            opt.read(blok);
            pp = execute_set(cat, name, opt);
            
            unsigned ix;
#ifdef BACKWARD_COMPATIBILITY
            // name changed to `property_number` on 10.12.2017
            if ( opt.set(ix, "property_number") || opt.set(ix, "property_index") )
#else
            if ( opt.set(ix, "property_number") )
#endif
            {
                if ( ix != pp->number() )
                    throw InvalidSyntax("Property number missmatch");
            }
        }
        else if ( cat == "simul" )
        {
            // adjust the name of the 'simul' property:
            if ( simul.prop->name() == "undefined" )
                simul.prop->rename(name);
            VLOG(" simul is named `" << name << "'\n");
            opt.read(blok);
            execute_change(name, opt);
        }
        else if ( do_change )
        {
            opt.read(blok);
            execute_change(name, opt);
        }
    }
    else
    {
        name = cat;
        //in this form, 'set' changes the value of an existing Property
#ifdef BACKWARD_COMPATIBILITY
        if ( spec )
        {
            //set CLASS:PARAMETER NAME { VALUE }
            is.get();
            para = Tokenizer::get_symbol(is);
            // Patch to accept 'set simul:display * {}':
            if ( cat == "simul" )
            {
                name = simul.prop->name();
                // skip the name
                Tokenizer::get_token(is);
            }
            else
                name = Tokenizer::get_token(is);
        }
        else
#endif
        {
            // set NAME PARAMETER { VALUE }
            para = Tokenizer::get_symbol(is);
        }
        
        // set NAME { PARAMETER = VALUE }
        blok = Tokenizer::get_block(is, '{');
        
        if ( blok.empty() )
            throw InvalidSyntax("syntax error");
        
        if ( do_change )
        {
            if ( para.empty() )
                opt.read(blok);
            else
                opt.define(para, blok);

            execute_change(name, opt);
        }
        else if ( para == "display" )
        {
            opt.define(para, blok);
            execute_change(name, opt, false);
        }
    }

    if ( pp && opt.warnings(std::cerr) )
        show_lines(is, ipos);
}

//------------------------------------------------------------------------------
/**
 Change the value of one (or more) parameters for property `NAME`.

     change NAME
     {
       PARAMETER = VALUE
       ...
     }
 
 Short syntax:
 
    change NAME { PARAMETER = VALUE }

The NAME should have been defined previously with the command `set`.
It is also possible to change all properties of a particular class:
 
     change all CLASS
     {
       PARAMETER = VALUE
       ...
     }

Examples:
 
    change system { viscosity = 0.5; }
    change system display { back_color=red; }
    change actin { rigidity = 1; }
    change microtubule { rigidity = 20 }
    change all fiber { confine = inside, 10; }
    change all fiber display { color = white; }

 */

void Parser::parse_change(std::istream& is)
{
    std::streampos ipos = is.tellg();
    bool change_all = false;
    
    std::string name = Tokenizer::get_symbol(is);
    std::string para;
    
    if ( name == "all" )
    {
        change_all = true;
        name = Tokenizer::get_symbol(is);
        if ( !simul.isPropertyClass(name) )
            throw InvalidSyntax("`"+name+"' is not a known class of object");
    }

#ifdef BACKWARD_COMPATIBILITY
    // Read formats anterior to 3.11.2017
    if ( is.peek() == ':' )
    {
        //change CLASS:PARAMETER NAME { VALUE }
        //change CLASS:PARAMETER * { VALUE }
        is.get();
        para = Tokenizer::get_symbol(is);
        std::string str = Tokenizer::get_token(is);
        if ( str == "*" )
            change_all = true;
        else
            name = str;
    }
    else
#endif
    {
        //change NAME PARAMETER { VALUE }
        para = Tokenizer::get_symbol(is);
        
#ifdef BACKWARD_COMPATIBILITY
        if ( is.peek() == '*' )
        {
            //change CLASS * { VALUE }
            is.get();
            change_all = true;
        }
        else if ( simul.findProperty(name, para) )
        {
            //change CLASS NAME { VALUE }
            name = para;
            para = "";
        }
#endif
    }

    //change NAME { VALUE }
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( blok.empty() )
        throw InvalidSyntax("syntax error");
    
    Glossary opt;
    if ( do_change )
    {
        if ( para.empty() )
            opt.read(blok);
        else
            opt.define(para, blok);
        
        if ( change_all )
            execute_change_all(name, opt);
        else
            execute_change(name, opt);
 
        if ( opt.warnings(std::cerr, ~0U) )
            show_lines(is, ipos);
    }
    else if ( para == "display" )
    {
        opt.define("display", blok);
        if ( change_all )
            execute_change_all(name, opt);
        else
            execute_change(name, opt, false);
    }
}

//------------------------------------------------------------------------------
/**
 The command `new` creates one or more objects with given specifications:
 
     new [MULTIPLICITY] NAME
     {
       position         = POSITION, [SPACE]
       direction        = DIRECTION
       orientation      = ROTATION, [ROTATION]
       mark             = INTEGER
       required         = INTEGER
     }
 
 The NAME should have been defined previously with the command `set`.\n

 The other parameters are:
 
 Parameter        | Type      | Default | Description                          |
 -----------------|-----------|---------|---------------------------------------
 MULTIPLICITY     | INTEGER   |   1     | number of objects.
 `position`       | POSITION  | random  | initial position of the object.
 `orientation`    | ROTATION  | random  | a rotation specified with respect to the object's center of gravity.
 `orientation[1]` | ROTATION  | none    | a rotation specified around the origin.
 `direction`      | DIRECTION | random  | specifies the direction of a fiber.
 `mark`           | INTEGER   |   0     | specifies a mark to be given to all objects created.
 `required`       | INTEGER   |   0     | minimum number of objects that should be created.
 

 Note that `position` only applies to movable objects, and `orientation` will have an effect only on objects that can be rotated. In addition, `position[1]` and `orientation[1]` are relevant only if `(MULTIPLICITY > 1)`, and do not apply to the first object.\n
 
 
 Short syntax:
 
     new [MULTIPLICITY] NAME ( POSITION )
 
 Shorter syntax:
 
     new [MULTIPLICITY] NAME
 
*/

void Parser::parse_new(std::istream& is)
{
    std::streampos ipos = is.tellg();
    Glossary opt;
    unsigned cnt = 1;
    Tokenizer::get_integer(is, cnt);
    std::string name = Tokenizer::get_symbol(is);
    
#ifdef BACKWARD_COMPATIBILITY
    // Read formats anterior to 3.11.2017
    if ( simul.isPropertyClass(name) )
    {
        std::string str = Tokenizer::get_symbol(is);
        if ( !str.empty() )
            name = str;
    }
#endif

    // Syntax sugar: () specify only position
    std::string blok = Tokenizer::get_block(is, '(');
    if ( blok.empty() )
    {
        blok = Tokenizer::get_block(is, '{');
        opt.read(blok);
    }
    else {
        opt.define("position", 0, blok);
    }
    
    if ( do_new  &&  cnt > 0 )
    {
        if ( opt.nb_keys() == 0 )
        {
            execute_new(name, cnt);
        }
        else
        {
            size_t nb_objects = simul.nbObjects();
            
            // syntax sugar, to specify the position of the Fiber ends
            if ( opt.has_key("position_ends") )
            {
                Vector A, B;
                if ( !opt.set(A, "position_ends") || !opt.set(B, "position_ends", 1) )
                    throw InvalidParameter("two vectors need to be defined by `position_ends'");
                opt.define("length",      0, (A-B).norm());
                opt.define("position",    0, (A+B)*0.5);
                opt.define("orientation", 0, (B-A).normalized());
            }

            // distribute objects regularly between two points:
            if ( opt.has_key("range") )
            {
                Vector A, B;
                if ( !opt.set(A, "range") || !opt.set(B, "range", 1) )
                    throw InvalidParameter("two vectors need to be defined by `range'");
                if ( opt.has_key("position") )
                    throw InvalidParameter("cannot specify `position' if `range' is defined");
                if ( cnt > 1 )
                {
                    Vector dAB = ( B - A ) / real(cnt-1);
                    for ( unsigned n = 0; n < cnt; ++n )
                    {
                        opt.define("position", 0, A + n * dAB);
                        execute_new(name, opt);
                    }
                }
                else
                {
                    opt.define("position", 0, 0.5*(A+B));
                    execute_new(name, opt);
                }
            }
            else
            {
                // place each object independently from the others:
                for ( unsigned n = 0; n < cnt; ++n )
                    execute_new(name, opt);
            }
            
            size_t required = 0;
            if ( opt.set(required, "required") )
            {
                size_t created = simul.nbObjects() - nb_objects;
                if ( created < required )
                {
                    std::cerr << "created  = " << created << '\n';
                    std::cerr << "required = " << required << '\n';
                    throw InvalidParameter("could not create enough `"+name+"'");
                }
            }
            
            if ( opt.has_key("display") )
                throw InvalidParameter("display parameters should be specified within `set'");
            
            if ( opt.warnings(std::cerr, ~0U) )
                show_lines(is, ipos);
        }
    }
}

//------------------------------------------------------------------------------
/**
 Delete objects:

     delete [MULTIPLICITY] NAME
     {
        mark      = INTEGER
        position  = [inside|outside], SPACE
        state     = [0|1], [0|1]
     }
 
 MULTIPLICITY is an integer, or the keyword 'all'.
 
 The parameters (mark, position, state) are all optional.
 All specified conditions must be fulfilled (this is a logical AND).
 The parameter `state` refers to bound/unbound state of Hands for Single and Couple,
 and to dynamic state for Fibers:
 - for Single, `state[0]` refers to the Hand: 0=free, 1=attached.
 - for Couple, `state[0]` refers to the first Hand: 0=free, 1=attached,
           and `state[1]` refers to the second Hand: 0=free, 1=attached.
 - for Fibers, `state[0]` refers to the Dynanic state of the PLUS end,
           and `state[1]` refers to the Dynanic state of the MINUS end.
 .
 
 To delete all objects of specified NAME:
 
     delete all NAME
 
 To delete at most COUNT objects of class NAME:
 
     delete COUNT NAME
 
 To delete all objects with a specified mark:
 
     delete all NAME
     {
        mark = INTEGER
     }
 
 To delete all objects within a Space:
 
     delete NAME
     {
        position = inside, SPACE
     }
 
 The SPACE must be the name of an existing Space.
 Only 'inside' and 'outside' are valid specifications.

 To delete all Couple called NAME that are not bound:
 
     delete all NAME { state = 0, 0; }
 
 To delete all Couple TYPE that are not crosslinking, use two calls:
 
     delete all NAME { state1 = 0; }
     delete all NAME { state2 = 0; }
 
 */

void Parser::parse_delete(std::istream& is)
{
    std::streampos ipos = is.tellg();
    unsigned cnt = 1;
    bool has_cnt = Tokenizer::get_integer(is, cnt);
    std::string name = Tokenizer::get_symbol(is);
#ifdef BACKWARD_COMPATIBILITY
    // Read formats anterior to 3.11.2017
    if ( simul.isPropertyClass(name) )
    {
        std::string str = Tokenizer::get_symbol(is);
        if ( !str.empty() )
            name = str;
    }
#endif
    if ( !has_cnt  &&  name == "all" )
    {
        cnt = ~0U; // this is very large
        name = Tokenizer::get_symbol(is);
    }
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_new )
    {
        Glossary opt(blok);
        execute_delete(name, opt, cnt);
        if ( opt.warnings(std::cerr) )
            show_lines(is, ipos);
    }
}


/**
 Mark objects:
 
     mark [MULTIPLICITY] NAME
     {
       mark       = INTEGER
       position   = POSITION
     }
 
 or
 
     mark all NAME
     {
         mark       = INTEGER
         position   = POSITION
     }

 NAME can be '*', and the parameter 'position' is optional.
 The syntax is the same as for command `delete`.
 */

void Parser::parse_mark(std::istream& is)
{
    std::streampos ipos = is.tellg();
    unsigned cnt = 0;
    bool has_cnt = Tokenizer::get_integer(is, cnt);
    std::string name = Tokenizer::get_symbol(is);
#ifdef BACKWARD_COMPATIBILITY
    // Read formats anterior to 3.11.2017
    if ( simul.isPropertyClass(name) )
    {
        std::string str = Tokenizer::get_symbol(is);
        if ( !str.empty() )
            name = str;
    }
#endif
    if ( !has_cnt  &&  name == "all" )
        name = Tokenizer::get_symbol(is);
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_new )
    {
        Glossary opt(blok);
        execute_mark(name, opt, cnt);
        if ( opt.warnings(std::cerr) )
            show_lines(is, ipos);
    }
}

//------------------------------------------------------------------------------
/**
 Cut all fibers that intersect a given plane.
 
     cut fiber NAME
     {
        plane = VECTOR, REAL
     }
 
     cut all fiber
     {
        plane = VECTOR, REAL
     }

 NAME can be '*' to cut all fibers.
 The plane is specified by a normal vector `n` (VECTOR) and a scalar `a` (REAL).
 The plane is defined by <em> n.pos + a = 0 </em>
 */

void Parser::parse_cut(std::istream& is)
{    
    std::streampos ipos = is.tellg();
    std::string name = Tokenizer::get_token(is);
 
    if ( name == "all" )
    {
        std::string str = Tokenizer::get_symbol(is);
        if ( str != "fiber" )
            throw InvalidParameter("only 'fiber' can be cut");
    }
    else if ( name == "fiber" )
    {
        name = Tokenizer::get_symbol(is);
    }
    
    std::string blok = Tokenizer::get_block(is, '{');

    if ( blok.empty() )
        throw InvalidSyntax("missing block after `cut'");
    
    if ( do_run )
    {
        Glossary opt(blok);
        execute_cut(name, opt);
        if ( opt.warnings(std::cerr) )
            show_lines(is, ipos);
    }
}


//------------------------------------------------------------------------------

/**
 @copydetails Interface::execute_run
 */
void Parser::parse_run(std::istream& is)
{
    std::streampos ipos = is.tellg();
    unsigned cnt = 1;
    bool has_cnt = Tokenizer::get_integer(is, cnt);
    std::string name = Tokenizer::get_symbol(is);
    
#ifdef BACKWARD_COMPATIBILITY
    // Read formats anterior to 3.11.2017
    if ( name == "simul" )
    {
        name = Tokenizer::get_symbol(is);
        if ( is.peek() == '*' )
        {
            is.get();
            name = simul.prop->name();
        }
    }
#endif
    if ( name.empty() )
        throw InvalidSyntax("missing simul name (use `run NB_STEPS all { }')");
    
    if ( name == "all" )
    {
        if ( Tokenizer::get_symbol(is) != "simul" )
            throw InvalidSyntax("expected `run all simul { }')");
        // There can only be one Simul object:
        name = simul.prop->name();
    }

    if ( name != simul.prop->name() )
        throw InvalidSyntax("unknown simul name `"+name+"'");

    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_run )
    {
        Glossary opt(blok);

        if ( !has_cnt )
        {
            // read `nb_steps' from the option block:
            if ( opt.set(cnt, "nb_steps") )
                opt.clear("nb_steps");
            else
            {
                // instead of `nb_steps', user can specify a duration in seconds:
                real span = 0;
                if ( opt.set(span, "duration") || opt.set(span, "time") )
                {
                    if ( span <= 0 )
                        throw InvalidParameter("duration must be >= 0'");
                    cnt = (unsigned)ceil(span/simul.prop->time_step);
                }
            }
        }
        else if ( opt.has_key("nb_steps") )
            throw InvalidSyntax("the number of steps was specified twice");
        
        if ( ! do_write )
            opt.clear("nb_frames");
        
        if ( opt.empty() )
            execute_run(cnt);
        else
            execute_run(cnt, opt);

        if ( opt.warnings(std::cerr) )
            show_lines(is, ipos);
    }
}

//------------------------------------------------------------------------------
/**
 Read and execute another config file.
 
     read FILE_NAME
     {
       required = BOOL
     }
 
 By default, `required = 1`, and execution will terminate if the file is not found.
 However, if `required=0`, the file will be executed if it is found, but execution
 will continue otherwise.
 
 \todo: able to specify `do_set` and `do_new` for command 'read'
*/

void Parser::parse_read(std::istream& is)
{
    std::streampos ipos = is.tellg();
    bool required = true;
    std::string file = Tokenizer::get_filename(is);
    
    if ( file.empty() )
        throw InvalidSyntax("missing/invalid file name after 'read'");
    
    std::string blok = Tokenizer::get_block(is, '{');
    if ( ! blok.empty() )
    {
        Glossary opt(blok);
        opt.set(required, "required");
        if ( opt.warnings(std::cerr) )
            show_lines(is, ipos);
    }
    
    if ( readConfig(file) )
    {
        if ( required )
            throw InvalidSyntax("could not open file `"+file+"'");
        else
            Cytosim::warn << "could not open file `" << file << "\n";
    }
}

//------------------------------------------------------------------------------
/**
 Import a simulation snapshot from a trajectory file
 
    import WHAT FILE_NAME
    {
        append = BOOL
        frame = INTEGER
    }
 
 The frame to be imported can be specified as an option: `frame=INTEGER`:
 
     import all my_file.cmo { frame = 10 }
 
 By default, this will replace the simulation state by the one loaded from file.
 To add the file objects to the simulation without deleting any of the current 
 object, you should specify `append = 1`:
 
     import all my_file.cmo { append = 1 }
 
 Finally instead of importing all the objects from the file, one can restrict
 the import to a desired class:
 
     import fiber my_file.cmo { append = 1 }
 
 Note that the simulation time will be changed to the one specified in the file,
 but this behavior can be changed by specifying the time:
 
     change system { time = 0 }
 
 ...assuming that the simul is called `system`.
 */

void Parser::parse_import(std::istream& is)
{
    std::streampos ipos = is.tellg();
    std::string what = Tokenizer::get_token(is);
    std::string file = Tokenizer::get_filename(is);
    
    if ( what.empty() )
        throw InvalidSyntax("missing class specification (use `import all FILENAME')");

    if ( file.empty() )
        throw InvalidSyntax("missing/invalid file name (use `import all FILENAME')");
    
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_new )
    {
        Glossary opt(blok);
        execute_import(file, what, opt);
        if ( opt.warnings(std::cerr) )
            show_lines(is, ipos);
    }
}


/**
 Export state to file. The general syntax is:
 
     export WHAT FILE_NAME
     {
       append = BOOL
       binary = BOOL
     }
 
 WHAT must be ``objects`` of ``properties``, and by default, both `binary` 
 and `append` are `true`. If `*` is specified instead of a file name,
 the current trajectory file will be used.
 
 Short syntax:
 
     export objects FILE_NAME
 
 Examples:
 
     export all sim_objects.cmo { append=0 }
     export properties properties.txt
 
 Attention: this command is disabled for `play`.
 */

void Parser::parse_export(std::istream& is)
{
    std::streampos ipos = is.tellg();
    std::string what = Tokenizer::get_token(is);
    std::string file = Tokenizer::get_filename(is);
    
    if ( what.empty() )
        throw InvalidSyntax("missing class specification (use `export all FILENAME')");
    
    if ( file.empty() )
        throw InvalidSyntax("missing/invalid file name (use `export all FILENAME')");
    

    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_write )
    {
        Glossary opt(blok);
        execute_export(file, what, opt);
        if ( opt.warnings(std::cerr) )
            show_lines(is, ipos);
    }
}


/**
 Export formatted data to file. The general syntax is:
 
     report WHAT FILE_NAME
     {
       append = BOOL
     }
 
 Short syntax:
 
     report WHAT FILE_NAME
 
 WHAT should be a valid argument to `report`:
 @copydetails Simul::report
 
 If `*` is specified instead of a file name, the report is sent to the standard output.
 
 Examples:
 
     report parameters parameters.cmo { append=0 }
     report fiber:length fibers_length.txt
 
 Note that this command is disabled for `play`.
 */

void Parser::parse_report(std::istream& is)
{
    std::streampos ipos = is.tellg();
    std::string what = Tokenizer::get_symbols(is);
    std::string file = Tokenizer::get_filename(is);

    if ( file.empty() )
        throw InvalidSyntax("missing/invalid file name after 'report'");
    
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_run && ( do_write || file == "*" ))
    {
        Glossary opt(blok);
        execute_report(file, what, opt);
        if ( opt.warnings(std::cerr) )
            show_lines(is, ipos);
    }
}

//------------------------------------------------------------------------------

/**
 Call custom function:
 
     call FUNCTION_NAME { OPTIONS }
 
 FUNCTION_NAME should be `equilibrate`, `custom0`, `custom1`, ... `custom9`.

 Note: The Simul::custom() functions need to be modified, to do something!
 */
void Parser::parse_call(std::istream& is)
{
    std::streampos ipos = is.tellg();
    std::string str = Tokenizer::get_symbol(is);
    
    if ( str.empty() )
        throw InvalidSyntax("missing function name after 'call'");
    
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_run )
    {
        Glossary opt(blok);
        execute_call(str, opt);
        if ( opt.warnings(std::cerr) )
            show_lines(is, ipos);
    }
}

//------------------------------------------------------------------------------
/**
 Repeat specified code.
 
     repeat INTEGER { CODE }
 
 */

void Parser::parse_repeat(std::istream& is)
{
    unsigned cnt = 1;
    
    if ( ! Tokenizer::get_integer(is, cnt) )
        throw InvalidSyntax("expected positive integer after 'repeat'");

    std::string code = Tokenizer::get_block(is, '{');
    
    for ( unsigned c = 0; c < cnt; ++c )
    {
        evaluate(code);
    }
}


/**
 Repeat specified code.
 
     for VAR=INTEGER:INTEGER { CODE }
 
 The two integers are the first and last iterations counters.
 Any occurence of VAR is replaced before the code is executed.
 
 Example:
 
     for CNT=1:10 {
       new 10 filament { length = CNT }
     }
 
 NOTE: This code is a hack, and can be improved vastly!
 */
void Parser::parse_for(std::istream& is)
{
    unsigned int start = 0, end = 1;
    
    std::string var = Tokenizer::get_symbol(is);
    
    std::string s = Tokenizer::get_token(is);
    if ( s != "=" )
        throw InvalidSyntax("missing '=' in command 'for'");
    
    if ( ! Tokenizer::get_integer(is, start) )
        throw InvalidSyntax("missing number after 'repeat'");

    s = Tokenizer::get_token(is);
    if ( s != ":" )
        throw InvalidSyntax("missing ':' in command 'for'");
    
    if ( ! Tokenizer::get_integer(is, end) )
        throw InvalidSyntax("missing number after 'repeat'");
    
    std::string code = Tokenizer::get_block(is, '{');
    
    for ( size_t c = start; c < end; ++c )
    {
        std::string sub = code;
        // substitute Variable name for this iteration:
        StreamFunc::find_and_replace(sub, var, std::to_string(c));
        //we use a fresh stream and Parser for each instance:
        Parser(*this).evaluate(sub);
        //hold();
    }
}

//------------------------------------------------------------------------------

/**
 Terminates execution
 
     end
 
 */
void Parser::parse_end(std::istream& is)
{
    if ( do_run )
        throw Exception("terminating program at command 'end'");

    /*
    std::string str = Tokenizer::get_symbol(is);
    
    if ( str == "if" )
    {
        str = Tokenizer::get_token(is);
        ABORT_NOW("the 'if' condition has not been implemented yet");
    }
     */
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Read and execute the next command to be found in the stream.
 Returns:
 - 0 if successfull
 - 1 if file is exhausted, or has error
 - 2 if 'end' was found.
 Thus parsing should be repeated while the return value is 0.
*/
int Parser::evaluate_one(std::istream& is)
{
    std::string tok = Tokenizer::get_token(is);
    
    if ( is.fail() )
        return 1;
    
    // empty lines
    if ( tok.length() < 1 || isspace(tok[0]) )
        return 0;
    
    // matlab-style comments (%{ })
    if ( tok[0] == '%' )
    {
        int c = is.get();
        if ( c == '{' )
            Tokenizer::get_block_text(is, 0, '}');
        else
            Tokenizer::get_line(is);
        return 0;
    }
    
    /*
     skip C-style comments:
     - single-line comment start with '/' and '/'
     - multi-line comments start with '/' and '*'
     */
    if ( tok[0] == '/' )
    {
        int c = is.get();
        if ( '/' == c )
            Tokenizer::get_line(is);
        else if ( '*' == c )
            Tokenizer::get_until(is, "*/");
        else
            throw InvalidSyntax("unexpected token after / `"+tok+"'");
        return 0;
    }
    
    //if (spos) StreamFunc::print_lines(std::clog, is, *spos, *spos);
    
    if ( tok == "set" )
        parse_set(is);
    else if ( tok == "change" )
        parse_change(is);
    else if ( tok == "new" || tok == "add" )
        parse_new(is);
    else if ( tok == "delete" )
        parse_delete(is);
    else if ( tok == "mark" )
        parse_mark(is);
    else if ( tok == "run" )
        parse_run(is);
    else if ( tok == "read" )
        parse_read(is);
#ifdef BACKWARD_COMPATIBILITY
    else if ( tok == "include" )
        parse_read(is);
#endif
    else if ( tok == "cut" )
        parse_cut(is);
    else if ( tok == "report" )
        parse_report(is);
#ifdef BACKWARD_COMPATIBILITY
    else if ( tok == "write" )
        parse_report(is);
#endif
    else if ( tok == "import" )
        parse_import(is);
    else if ( tok == "export" )
        parse_export(is);
    else if ( tok == "call" )
        parse_call(is);
    else if ( tok == "repeat" )
        parse_repeat(is);
    else if ( tok == "skip" )
        Tokenizer::get_block(is, '{');
    else if ( tok == "for" )
        parse_for(is);
    else if ( tok == "restart" )
    {
        static unsigned cnt = 0;
        unsigned num = 1;
        Tokenizer::get_integer(is, num);
        if ( do_run && cnt++ < num )
        {
            // reset simulation and rewind config file
            simul.erase();
            is.clear();
            is.seekg(0);
            return 0;
        }
    }
    else if ( tok == "end" )
        return 2; //parse_end(is);
    else if ( tok == ";" )
        return 0;
    else if ( tok == "dump" )
    {
        if ( do_write && do_run )
            simul.dump();
    }
    else if ( tok == "dump_system" )
    {
        if ( do_write && do_run )
            simul.dump_system();
    }
    else {
        throw InvalidSyntax("unexpected command `"+tok+"'");
    }
    return 0;
}


void Parser::evaluate(std::istream& is)
{
    std::streampos ipos(0);
    try {
        while ( is.good() )
        {
            ipos = is.tellg();
            if ( evaluate_one(is) )
                return;
        }
    }
    catch( Exception & e )
    {
        //e << ", " + msg + ":\n";
        e << ":\n";
        e << StreamFunc::get_lines(is, ipos, is.tellg());
        throw;
    }
}


void Parser::evaluate(std::string const& code)
{
    std::istringstream is(code);
    evaluate(is);
}


int Parser::readConfig(std::string const& filename)
{
    std::ifstream is(filename.c_str(), std::ifstream::in);
    if ( is.good() )
    {
        VLOG(" Parsing `" << filename << "' set " << do_set << " change " << do_change);
        VLOG(" new " << do_new << " run " << do_run << " write " << do_write << " )\n");
        evaluate(is);
        return 0;
    }
    return 1;
}


int Parser::readConfig()
{
    return readConfig(simul.prop->config_file);
}

