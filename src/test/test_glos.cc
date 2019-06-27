// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "ansi_colors.h"
#include "exceptions.h"
#include "glossary.h"

int i = 0;
float f = 0.0;
std::string s;


int main(int argc, char* argv[])
{
    Glossary glos;
    
    try {
        // read command line arguments:
        glos.read_strings(argc-1, argv+1);
        
        // read file if provided on command line:
        // the file name is recognized by its extension
        std::string str;
        if ( glos.set(str, ".cym") )
            glos.read_file(str);
        
        // print content of Glossary:
        printf("%lu keys:\n", glos.nb_keys());
        glos.write(std::cout, "    > ");
        
        // extract values from Glossary:
        if ( glos.set(i, "integer") )  printf("integer : %i\n", i);
        if ( glos.set(f, "float") )    printf("float : %f\n", f);
        if ( glos.set(s, "string") )   printf("string : %s\n", s.c_str());
        
        return EXIT_SUCCESS;
    }
    catch ( Exception& e )
    {
        print_magenta(std::cout, "Error : "+e.message());
        return EXIT_FAILURE;
    }
}
