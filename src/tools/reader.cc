// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/*
 This is mostly a test for the class defined in frame_reader.h
 but it can be used to navigate from frame to frame in a object-file
*/

#include <cstring>
#include <cctype>
#include <cstdlib>

#include "glossary.h"
#include "messages.h"
#include "iowrapper.h"
#include "frame_reader.h"
#include "simul.h"
#include "parser.h"

//------------------------------------------------------------------------------

void help(std::ostream& os)
{
    printf("Cytosim-reader %iD, file version %i\n", DIM, Simul::currentFormatID);
    os << "\n";
    os << "Syntax:  reader [OPTIONS] INPUT_FILE_NAME output=FILE_NAME\n";
    os << "\n";
    os << "OPTIONS:\n";
    os << "     help       display this message\n";
    os << "     binary=0   write text coordinates in `file_out'\n";
    os << "     binary=1   write binary coordinates in `file_out'\n";
    os << "     verbose=?  set the verbose level\n";
    os << "\n";
}

void instructions(std::ostream& os = std::cout)
{
    os << "Commands understood at prompt:\n";
    os << "  'q'      quit\n";
    os << "  'n'      read next frame\n";
    os << "  'w'      write frame\n";
    os << "  'c'      clear buffer without changing positions\n";
    os << "  'r'      rewind\n";
    os << "  'e'      erase state\n";
    os << " INTEGER   read specified frame if possible\n";
}


//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    int frame;
    Simul simul;
    char user[1024] = "\0";
    FrameReader reader;
    Glossary arg;
    
    arg.read_strings(argc-1, argv+1);
    
    if ( arg.use_key("help") )
    {
        help(std::cout);
        instructions();
        return EXIT_SUCCESS;
    }
    
    std::string input = TRAJECTORY;
    std::string output = "output.cmo";
    
    arg.set(output, "output");
    arg.set(input, "input") || arg.set(input, ".cmo");

    bool binary = true;
    arg.set(binary, "binary");
    
    try {
        simul.loadProperties();
        reader.openFile(input);
    }
    catch( Exception & e )
    {
        std::clog << "Aborted: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    if ( !reader.good() )
    {
        printf("File could not be oppened\n");
        return EXIT_FAILURE;
    }
    
    
    printf("Reader: enter (h) for help\n");
    while ( true )
    {
        if ( reader.hasFrame() )
        {
            printf("Reader: Frame %li in buffer\n", reader.currentFrame());
            simul.reportInventory(std::cout);
        }
        else
            printf("Reader: Empty buffer\n");
        
        printf("Reader: ");
        fgets(user, sizeof(user), stdin);
        
        if ( isdigit( user[0] ))
        {
            if ( 1 == sscanf(user, "%i", &frame ) )
            {
                try {
                    if ( 0 != reader.loadFrame(simul, frame) )
                        printf("Reader: frame not found: ");
                }
                catch( Exception & e ) {
                    printf("Reader: exception in `read` %i: %s\n", frame, e.c_str());
                }
            }
        }
        else
        {
            switch( user[0] )
            {
                case '\n':
                case 'n':
                    try {
                        int err = reader.loadNextFrame(simul);
                        if ( err ) printf("Reader error with `next`: %i\n", err);
                    }
                    catch( Exception & e ) {
                        printf("Reader: exception in `next`: %s\n", e.c_str());
                    }
                    break;
                    
                case 'w':
                    simul.writeObjects(output, true, binary);
                    break;
                    
                case 'e':
                    simul.erase();
                    break;
                    
                case 'b':
                    binary = !binary;
                    printf("Reader: binary = %i\n", binary);
                    break;
                    
                case 'c':
                    reader.clearPositions();
                    break;

                case 'r':
                    reader.rewind();
                    break;
                    
                case 'q': case 'Q': case 27:
                    return 0;
                    
                default:
                    instructions();
                    break;
            }
        }
    }
}
