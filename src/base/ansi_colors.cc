// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on 03/09/2014.


#include "ansi_colors.h"
#include <cstdio>
#include <cstdlib>


#if ( 0 )


/* Return the number of colors that the terminal supports */
int nb_colors_supported()
{
    unsigned long n = 0;
    FILE * fp = popen("tput colors 2> /dev/null", "r");
    if ( fp ) {
        char str[32] = { 0 };
        if ( fgets(str, sizeof(str), fp) )
            n = strtoul(str, nullptr, 10);
        printf("tput: %s\n", str);
        pclose(fp);
    }
    //printf("%i colors\n", n);
    return (int)n;
}


void print_red(std::ostream& os, std::string const& str)
{
    if ( nb_colors_supported() > 7 )
        os << KBLDRED << str << KNRM;
    else
        os << str;
}

void print_green(std::ostream& os, std::string const& str)
{
    if ( nb_colors_supported() > 7 )
        os << KBLDGRN << str << KNRM;
    else
        os << str;
}

void print_yellow(std::ostream& os, std::string const& str)
{
    if ( nb_colors_supported() > 7 )
        os << KBLDYEL << str << KNRM;
    else
        os << str;
}

void print_blue(std::ostream& os, std::string const& str)
{
    if ( nb_colors_supported() > 7 )
        os << KBLDBLU << str << KNRM;
    else
        os << str;
}

void print_magenta(std::ostream& os, std::string const& str)
{
    if ( nb_colors_supported() > 7 )
        os << KBLDMAG << str << KNRM;
    else
        os << str;
}

void print_cyan(std::ostream& os, std::string const& str)
{
    if ( nb_colors_supported() > 7 )
        os << KBLDCYN << str << KNRM;
    else
        os << str;
}

#else


void print_red(std::ostream& os, std::string const& str)
{
    os << str;
}

void print_green(std::ostream& os, std::string const& str)
{
    os << str;
}

void print_yellow(std::ostream& os, std::string const& str)
{
    os << str;
}

void print_blue(std::ostream& os, std::string const& str)
{
    os << str;
}

void print_magenta(std::ostream& os, std::string const& str)
{
    os << str;
}

void print_cyan(std::ostream& os, std::string const& str)
{
    os << str;
}

#endif

