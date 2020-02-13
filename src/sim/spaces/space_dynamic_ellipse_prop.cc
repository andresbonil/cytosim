// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <cmath>
#include "space_dynamic_ellipse_prop.h"
#include "space_dynamic_ellipse.h"
#include "space_prop.h"
#include "property_list.h"
#include "glossary.h"


//------------------------------------------------------------------------------


Space * SpaceDynamicEllipseProp::newSpace() const
{
	return new SpaceDynamicEllipse(this);
}


Space * SpaceDynamicEllipseProp::newSpace(Glossary& opt) const
{
    Space * spc = newSpace();
    
    if ( spc )
    {
        // normal way to set the size:
        spc->resize(opt);
    }
    return spc;
}


void SpaceDynamicEllipseProp::clear()
{
    tension = 0 ;
    volume  = 0  ;
	SpaceProp::clear();
}

void SpaceDynamicEllipseProp::read(Glossary& glos)
{
    SpaceProp::read(glos);
	glos.set(tension, "tension");

}



void SpaceDynamicEllipseProp::complete(Simul const& sim) 
{
	SpaceProp::complete(sim);
	
	if	(tension < 0)
		throw InvalidParameter("tension must be positive");
}


void SpaceDynamicEllipseProp::write_values(std::ostream& os) const
{
    SpaceProp::write_values(os);
	write_value(os, "tension", tension);
}