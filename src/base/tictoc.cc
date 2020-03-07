// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on 27/9/2008.


#include "tictoc.h"
#include <cstdio>
#include <ctime>
#include <cstring>
#include <sys/time.h>

#pragma mark Wall time

/**
 This get current time from the C-library functions time() and ctime_r()
 */
void TicToc::get_date(char * buf, size_t buf_size)
{
    if ( buf_size > 25 )
    {
        time_t now = time(nullptr);
        //asctime_r(localtime(&now), buf);
#if ( 0 )
        ctime_r(&now, buf);
#else
        strncpy(buf, ctime(&now), buf_size);
#endif
        // remove new line:
        buf[24] = 0;
    }
    else if ( buf_size > 0 )
        buf[0] = 0;
    // terminate string:
    if ( buf_size > 0 )
        buf[buf_size-1] = 0;
}


void TicToc::get_date(char * buf, size_t buf_size, bool no_year)
{
    if ( buf_size > 25 )
    {
        TicToc::get_date(buf, buf_size);
        // remove year:
        if ( no_year )
            buf[19] = 0;
    }
}


char const* TicToc::date()
{
    static char buf[32];
    get_date(buf, sizeof(buf));
    return buf;
}


int  TicToc::days_since_2000()
{
    time_t now = time(nullptr);
    tm * loc = localtime(&now);
    return loc->tm_year - 2000 + loc->tm_yday;
}


time_t TicToc::seconds_since_1970()
{
    return time(nullptr);
}

time_t  TicToc::seconds_since_2000()
{
    return time(nullptr) - 946684800;
}


int  TicToc::year()
{
    time_t now = time(nullptr);
    tm * loc = localtime(&now);
    return loc->tm_year;
}


int  TicToc::day_of_the_year()
{
    time_t now = time(nullptr);
    tm * loc = localtime(&now);
    return loc->tm_yday;
}


int  TicToc::hours_today()
{
    time_t now = time(nullptr);
    tm * loc = localtime(&now);
    return loc->tm_hour;
}


double TicToc::seconds_today()
{
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return (double)tv.tv_sec + 1e-6 * tv.tv_usec;
}


double TicToc::centiseconds()
{
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return 100 * tv.tv_sec + tv.tv_usec / 10000.0;
}

double TicToc::microseconds()
{
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return tv.tv_usec;
}

#pragma mark CPU time


/**
 This calls the C-library function clock()
 */
double TicToc::processor_time(clock_t& clk)
{
    clock_t now = clock();
    double sec = double( now - clk ) / CLOCKS_PER_SEC;
    clk = now;
    return sec;
}


#if ( 1 )

// using real time
struct timeval tic_t;

const int nbins = 8;
double bins[nbins] = { 0 };

void TicToc::tic()
{
    gettimeofday(&tic_t, nullptr);
}

double TicToc::toc()
{
    timeval tv;
    gettimeofday(&tv, nullptr);
    return (tv.tv_sec-tic_t.tv_sec) + 1e-6*(tv.tv_usec-tic_t.tv_usec);
}

void TicToc::toc(int n)
{
    if ( n < nbins )
    {
        timeval tv;
        gettimeofday(&tv, nullptr);
        bins[n] += tv.tv_sec-tic_t.tv_sec + 1e-6*(tv.tv_usec-tic_t.tv_usec);
        tic_t = tv;
    }
}

void TicToc::report(const char * msg)
{
    for ( int n = 0; n < nbins; ++n )
    {
        printf("%s %i : %6.3fs\n", msg, n, bins[n]);
        bins[n] = 0;
    }
}

#else

#include <ctime>
clock_t tic_t;

// using the CPU time
void TicToc::tic()
{
    tic_t = clock();
}

double TicToc::toc()
{
    return (double)( clock() - tic_t ) / CLOCKS_PER_SEC;
}

#endif


void TicToc::toc(const char * msg)
{
    printf("%s : %6.6fs\n", msg, toc());
}


void TicToc::toc(const char * msg, const char * end)
{
    if ( end )
        printf("%s : %6.3fs %s", msg, toc(), end);
    else
        printf("%s : %6.3fs ", msg, toc());
}


