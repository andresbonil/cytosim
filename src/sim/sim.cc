// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "simul.h"
#include "parser.h"
#include "messages.h"
#include "glossary.h"
#include "exceptions.h"
#include "backtrace.h"
#include "splash.h"
#include "tictoc.h"
#include <csignal>
#include "unistd.h"


void help(std::ostream& os)
{
    os << "sim [OPTIONS] [FILE]\n";
    os << "  FILE    run specified config file (FILE must end with `.cym')\n";
    os << "  *       print messages to terminal (and not `messages.cmo')\n";
    os << "  info    print build options\n";
    os << "  help    print this message\n";
}


void handle_signal(int sig)
{
    /*
     A signal handler is restricted to call only async-signal-safe-functions
     practically speaking, most syscalls(2) only
     */
    char str[128] = { 0 };
    strncpy(str, "Cytosim received signal   \n", 128);
    str[25] = '0' +  sig     % 10;
    str[24] = '0' + (sig/10) % 10;
    write(STDERR_FILENO, str, 28);
    _exit(sig);
}


void handle_interrupt(int sig)
{
    Cytosim::out << "killed " << sig << "\n" << std::endl;
    _exit(sig);
}

//------------------------------------------------------------------------------
//=================================  MAIN  =====================================
//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    // register callback to catch interrupting signals:
    if ( signal(SIGINT, handle_interrupt) )
        std::cerr << "Could not register SIGINT handler\n";
    if ( signal(SIGTERM, handle_interrupt) )
        std::cerr << "Could not register SIGTERM handler\n";
    if ( signal(SIGSEGV, handle_signal) )
        std::cerr << "Could not register SIGSEGV handler\n";
    if ( signal(SIGILL,  handle_signal) )
        std::cerr << "Could not register SIGILL handler\n";
    if ( signal(SIGABRT, handle_signal) )
        std::cerr << "Could not register SIGABRT handler\n";

    Glossary arg;

    //parse the command line:
    arg.read_strings(argc-1, argv+1);
    
    if ( arg.use_key("help") || arg.use_key("--help") )
    {
        splash(std::cout);
        help(std::cout);
        return EXIT_SUCCESS;
    }

    if ( arg.use_key("info") || arg.use_key("--version")  )
    {
        print_version(std::cout);
        std::cout << "    DIM = " << DIM << '\n';
        return EXIT_SUCCESS;
    }
    
    if ( ! arg.use_key("*") )
    {
        Cytosim::out.open("messages.cmo");
        Cytosim::log.redirect(Cytosim::out);
        Cytosim::warn.redirect(Cytosim::out);
    }
    
#ifdef CODE_VERSION
    Cytosim::out << "CYTOSIM PI version " << CODE_VERSION << "\n";
#else
    Cytosim::out << "CYTOSIM PI\n";
#endif

    Simul simul;
    try {
        simul.initialize(arg);
    }
    catch( Exception & e ) {
        std::cerr << "Error: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
    catch(...) {
        std::cerr << "Error: an unknown exception occured during initialization\n";
        return EXIT_FAILURE;
    }
    
    arg.warnings(std::cerr);
    time_t sec = TicToc::seconds_since_1970();
    
    try {
        if ( Parser(simul, 1, 1, 1, 1, 1).readConfig() )
            std::cerr << "You must specify a config file\n";
    }
    catch( Exception & e ) {
        std::cerr << "\nError: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
    catch(...) {
        std::cerr << "\nAn unknown exception occured\n";
        return EXIT_FAILURE;
    }
    
    Cytosim::out << "% " << TicToc::date() << "\n";
    sec = TicToc::seconds_since_1970() - sec;
    Cytosim::out << "end  " << sec << " s ( " << ( sec / 60 ) / 60.0 << " h )\n";
    Cytosim::out.close();
    return EXIT_SUCCESS;
}