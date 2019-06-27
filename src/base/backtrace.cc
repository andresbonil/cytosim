// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on 24/04/2010.


#include "backtrace.h"

// You can disable backtrace by changing the line below into '#if ( 0 )'
#ifdef __GNUC__


#include <execinfo.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

/**
 * print the current call-stack using functions from the GNU C Library:
 * - backtrace()
 * - backtrace_symbols()
 * .
 * provided by <execinfo.h>
 */
void print_backtrace(int out)
{
    void* buffer[128];
    int size = backtrace(buffer, 128);
#if ( 1 )
    char** strs = backtrace_symbols(buffer, size);

    write(out, "Cytosim execution stack:\n", 25);
    for ( int ii = 0; ii < size; ++ii )
    {
        write(out, strs[ii], strlen(strs[ii]));
        write(out, "\n", 1);
    }
    free(strs);
#else
    backtrace_symbols_fd(buffer, size, out);
#endif
}

#else

void print_backtrace(int out)
{
    write(out, "Execution stack information unavailable\n", 40);
}

#endif

