// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on 03/09/2014.


#include <iostream>

#define KNRM  "\x1B[0m"
#define KBLD  "\x1B[1m"
#define KUND  "\x1B[4m"
#define KREV  "\x1B[7m"

#define KBLK  "\x1B[30m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

#define KBLDRED  "\x1B[1m\x1B[31m"
#define KBLDGRN  "\x1B[1m\x1B[32m"
#define KBLDYEL  "\x1B[1m\x1B[33m"
#define KBLDBLU  "\x1B[1m\x1B[34m"
#define KBLDMAG  "\x1B[1m\x1B[35m"
#define KBLDCYN  "\x1B[1m\x1B[36m"
#define KBLDWHT  "\x1B[1m\x1B[37m"


/// Returns the number of colors that the terminal supports
int nb_colors_supported();

/// print text in red
void print_red(std::ostream&, std::string const&);

/// print text in green
void print_green(std::ostream&, std::string const&);

/// print text in yellow
void print_yellow(std::ostream&, std::string const&);

/// print text in blue
void print_blue(std::ostream&, std::string const&);

/// print text in magenta
void print_magenta(std::ostream&, std::string const&);

/// print text in cyan
void print_cyan(std::ostream&, std::string const&);


