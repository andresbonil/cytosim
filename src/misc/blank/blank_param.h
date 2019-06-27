// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on 6/27/07.

#ifndef BLANK_PARAM_H
#define BLANK_PARAM_H

#include "real.h"
#include <string>
class Glossary;

///The list of parameters for Blank
class BlankParam
{
public:
    
    //--------------------------------------------------------------------------

    unsigned int max;
    real         diffusion, diffusion_dt;
    
    real         time_step;
    real         time;
    
    unsigned int delay;
    unsigned int repeat;
    std::string  config;     //the parameter file
    

    //--------------------------------------------------------------------------
    
    
    /// constructor
    BlankParam()  { clear(); }
    
    /// set default values
    void clear();
    
    /// read
    void read(Glossary& glos);

    /// write all values
    void write(std::ostream&) const;
};

#endif
