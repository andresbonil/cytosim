// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on 27/9/2008.

#ifndef TICTOC_H
#define TICTOC_H

#include <ctime>

/// import __rdtsc()
#if defined(__APPLE__) || defined(__linux)
#  include <x86intrin.h>
#else
#  include <intrin.h>
#endif

/// A set of functions related to time
/**
 Functions to get wall-time, and processor-time,
 calling the C-standard library.
 */
namespace TicToc
{    
    
    /// set current date in short format, `buf` should be 26 character long or more
    void    get_date(char * buf, size_t buf_size);
    
    /// set current date in short format, `buf` should be 26 character long or more
    void    get_date(char * buf, size_t buf_size, bool no_year);

    /// using a local char[] to call get_date()
    char const* date();

    
    /// approximately the number of days since Jan 1 2000
    int     days_since_2000();
    
    /// number of seconds since Jan 1 1970, 0h00
    time_t  seconds_since_1970();
    
    /// number of seconds since Jan 1 2000, 0h00
    time_t  seconds_since_2000();

    /// year
    int     year();

    /// day of the year (0-365)
    int     day_of_the_year();
    
    /// hour of the day (0-23)
    int     hours_today();

    /// number of seconds since midnight
    double  seconds_today();
    
    /// number of centiseconds since midnight
    double  centiseconds();
 
    /// number of microseconds
    double  microseconds();
    
    /// return CPU time in seconds and update `clk` to current time
    double  processor_time(clock_t&);
    
    /// call to start timer
    void    tic();
    
    /// return number of seconds elapsed since last call to `tic()`
    double  toc();

    /// print the time elapsed since 'tic', as "msg: time"
    void    toc(const char * msg);
    
    /// print the time elapsed since 'tic', as "msg: time end"
    void    toc(const char * msg, const char * end);
    
    /// add elapsed time since last call to `tic()` to register 'n'
    void    toc(int);
    
    /// report values of all registers
    void    report(const char * msg);

}

#endif
