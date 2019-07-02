// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


/**
 Calculate the grid size automatically for dynamic Fiber,
 Bead, Sphere and Solid.
 
 This function can be used to set SimulProp::steric_max_range.
 
 We assume that Fiber::adjustSegmentation() is used, ensuring that
 ( actual segmentation ) < ( 4/3 * FiberProp::segmentation ).
 */
real Simul::estimateStericRange() const
{
    real ran = 0;
    real len = 0;
    
    // check all FiberProp with enabled steric:
    for ( Property * i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp*>(i);
        if ( fp->steric )
        {
            // The maximum length of a segment is 4/3 * segmentation
            len = std::max(len, (real)(1.4) * fp->segmentation);
            
            // check extended range of interaction
            ran = std::max(ran, fp->steric_radius + fp->steric_range);
        }
    }
    
    // verify against the actual segmentations of the Fibers:
    for ( Fiber const* fib=fibers.first(); fib; fib=fib->next() )
    {
        if ( fib->prop->steric )
            len = std::max(len, fib->segmentation());
    }

    /*
     The interaction can be aligned with the fiber, and we must add the distances:
     2 * range if two fibers of radius 'range' interact.
     + 2 * ( len / 2 ) since len/2 is the distance between the center of the segment
     and its most distal point.
     */
    ran =  len + 2*ran;
    
    
    for ( Sphere const* sp=spheres.first(); sp; sp=sp->next() )
    {
        if ( sp->prop->steric )
            ran = std::max(ran, 2 * sp->radius() + sp->prop->steric_range);
    }
    
    for ( Bead const* bd=beads.first(); bd; bd=bd->next() )
    {
        if ( bd->prop->steric )
            ran = std::max(ran, 2 * bd->radius() + bd->prop->steric_range);
    }
    
    for ( Solid const* so=solids.first(); so; so=so->next() )
    {
        if ( so->prop->steric )
        {
            for ( unsigned p = 0; p < so->nbPoints(); ++p )
                ran = std::max(ran, 2 * so->radius(p) + so->prop->steric_range);
        }
    }
    
    if ( ran < REAL_EPSILON )
        PRINT_ONCE("Warning: could not estimate simul:steric_max_range automatically!\n");
    
    return ran;
}


void Simul::setStericGrid(Space const* spc) const
{
    assert_true(spc);
    real& range = prop->steric_max_range;
    
    if ( range <= 0 ) 
    {
        range = estimateStericRange();
        //Cytosim::log("auto setting simul:steric_max_range=%.3f\n", range);
    }
    
    if ( range <= 0 )
        throw InvalidParameter("simul:steric_max_range must be defined");

    const size_t sup = 1 << 17;
    while ( pointGrid.setGrid(spc, range) > sup )
    {
        //std::clog << "increasing simul:steric_max_range\n";
        range *= 2;
    }
    pointGrid.createCells();
}


/**
 The prop->steric of each object is a bit-field that
 specify one or more 'pane' where the object is present.
 The different panes are then treated consecutively and independently, 
 and only objects in the same pane may interact.
 
     for ( int pane=1; pane<=2 && pane<=prop->steric; ++pane )
     {
         if ( obj->prop->steric & pane )
         ...
     }
 
 With this mechanism, the user can flexibly configure which objects
 may see each other and thus control the steric interactions.
 
 At present, we only support 1 pane (Simul property steric).
 This can be extended if necessary, but the steric_stiffness[]
 properties should be extended as well.
 */
void Simul::setStericInteractions(Meca& meca) const
{
    if ( !pointGrid.hasGrid() )
    {
        if (!spaces.master())
            return;
        setStericGrid(spaces.master());
    }

    // clear grid
    pointGrid.clear();
    
    // distribute Fiber-points on the grid
    for ( Fiber* fib=fibers.first(); fib; fib=fib->next() )
    {
        if ( fib->prop->steric )
        {
            const real rad = fib->prop->steric_radius;        // equilibrium radius
            const real ran = rad + fib->prop->steric_range;   // extended range of interaction
        
            // include segments, in the cell associated with their center
            for ( unsigned r = 0; r < fib->nbSegments(); ++r )
#if ( NB_STERIC_PANES == 1 )
                pointGrid.add(FiberSegment(fib, r), rad, ran);
#else
                pointGrid.add(fib->prop->steric, FiberSegment(fib, r), rad, ran);
#endif
        }
    }
    
    // include Spheres
    for ( Sphere* sp=spheres.first(); sp; sp=sp->next() )
    {
        if ( sp->prop->steric )
#if ( NB_STERIC_PANES == 1 )
            pointGrid.add(Mecapoint(sp, 0), sp->radius(), sp->radius()+sp->prop->steric_range);
#else
            pointGrid.add(sp->prop->steric, Mecapoint(sp, 0), sp->radius(), sp->radius()+sp->prop->steric_range);
#endif
    }
    
    // include Beads
    for ( Bead* bd=beads.first(); bd; bd=bd->next() )
    {
        if ( bd->prop->steric )
#if ( NB_STERIC_PANES == 1 )
            pointGrid.add(Mecapoint(bd, 0), bd->radius(), bd->radius()+bd->prop->steric_range);
#else
            pointGrid.add(bd->prop->steric, Mecapoint(bd, 0), bd->radius(), bd->radius()+bd->prop->steric_range);
#endif
    }
        
    // include Points that have a radius from Solids
    for ( Solid* so=solids.first(); so; so=so->next() )
    {
        if ( so->prop->steric )
        {
            for ( unsigned i = 0; i < so->nbPoints(); ++i )
            {
                if ( so->radius(i) > REAL_EPSILON )
#if ( NB_STERIC_PANES == 1 )
                    pointGrid.add(Mecapoint(so, i), so->radius(i), so->radius(i)+so->prop->steric_range);
#else
                    pointGrid.add(so->prop->steric, Mecapoint(so, i), so->radius(i), so->radius(i)+so->prop->steric_range);
#endif
            }
        }
    }
    
    /// create parameters
    PointGridParam pam(prop->steric_stiffness_push[0], prop->steric_stiffness_pull[0]);
    
#if ( NB_STERIC_PANES == 1 )
    
    pointGrid.setInteractions(meca, pam);

#elif ( NB_STERIC_PANES == 2 )
    
    // add steric interactions inside pane 1:
    pointGrid.setInteractions(meca, pam, 1);
    // add steric interactions between panes 1 and 2:
    pointGrid.setInteractions(meca, pam, 1, 2);
    //pointGrid.setInteractions(meca, pam, 2, 1);

#else
    
    // add steric interactions between different panes:
    for ( unsigned p = 1; p <= NB_STERIC_PANES; ++p )
        pointGrid.setInteractions(meca, pam, p);

#endif
}


//------------------------------------------------------------------------------
/**
 This will:
 - Register all Mecables in the Meca: Fiber Solid Bead and Sphere
 - call setInteractions() for all objects in the system,
 - call setStericInteractions() if prop->steric is true.
 .
 */
void Simul::setInteractions(Meca & meca) const
{
    // prepare the meca, and register Mecables
    meca.clear();
    
    for ( Fiber  * f= fibers.first(); f ; f=f->next() )
        meca.add(f);
    for ( Solid  * s= solids.first(); s ; s=s->next() )
        meca.add(s);
    for ( Sphere * o=spheres.first(); o ; o=o->next() )
        meca.add(o);
    for ( Bead   * b=  beads.first(); b ; b=b->next() )
        meca.add(b);
    
    meca.prepare(prop);
    
    // add interactions for all objects:
    
    for ( Space * s=spaces.first(); s; s=s->next() )
        s->setInteractions(meca, fibers);
    
    for ( Fiber * f=fibers.first(); f ; f=f->next() )
        f->setInteractions(meca);
    
    for ( Solid * s=solids.first(); s ; s=s->next() )
        s->setInteractions(meca);
    
    for ( Sphere * o=spheres.first(); o ; o=o->next() )
        o->setInteractions(meca);
    
    for ( Bead * b=beads.first(); b ; b=b->next() )
        b->setInteractions(meca);

    for ( Single * i=singles.firstA(); i ; i=i->next() )
        i->setInteractions(meca);

    for ( Couple * c=couples.firstAA(); c ; c=c->next() )
        c->setInteractions(meca);
    
    for ( Organizer * a = organizers.first(); a; a=a->next() )
        a->setInteractions(meca);

    //for ( Event * e = events.first(); e; e=e->next() )
    //    e->setInteractions(meca);

    // add steric interactions
    if ( prop->steric )
        setStericInteractions(meca);
    
    
    // ALL THE FORCES BELOW WERE DONE FOR TESTING PURPOSES:
#if ( 0 )
    PRINT_ONCE("AD-HOC CALIBRATED FORCE ENABLED\n");
    // add calibrated forces, for testing rotation
    for ( Fiber * fib = fibers.first(); fib; fib = fib->next() )
        meca.addTorqueClamp(fib->interpolateCenter(), Vector(0,1,0), 1);
#endif
#if ( 0 )
    PRINT_ONCE("AD-HOC CALIBRATED FORCE ENABLED\n");
    // add calibrated force to test rotation of spheres:
    Vector force(0,1,0);
    for ( Sphere * sph = spheres.first(); sph; sph = sph->next() )
    {
        meca.addForce(Mecapoint(sph, 1), -force);
        meca.addForce(Mecapoint(sph, 2), +force);
    }
#endif
}


/// solve the system
void Simul::solve()
{
    //auto rdtsc = __rdtsc();
    setInteractions(sMeca);
    //printf("     ::set      %16llu\n", (__rdtsc()-rdtsc)>>5); rdtsc = __rdtsc();
    sMeca.solve(prop, prop->precondition);
    //printf("     ::solve    %16llu\n", (__rdtsc()-rdtsc)>>5); rdtsc = __rdtsc();
    sMeca.apply();
    //printf("     ::apply    %16llu\n", (__rdtsc()-rdtsc)>>5);
}


/*
 Solve the system, and automatically select the fastest preconditionning method
 */
void Simul::solve_auto()
{
    setInteractions(sMeca);
    
    // solve the system, recording time:
    long cpu = TicToc::centiseconds();
    sMeca.solve(prop, precondMethod);
    cpu = TicToc::centiseconds() - cpu;
    
    sMeca.apply();

    // Automatic selection of preconditionning method:
    const unsigned N_TEST = 6;
    const unsigned PERIOD = 32;
    
    //automatically select the preconditionning mode:
    //by trying each methods N_STEP steps, adding CPU time and use fastest.
    if ( precondCounter++ < N_TEST )
    {
        precondCPU[precondMethod] += cpu;
        
        //std::clog << " precond "<<precondMethod<<" iter " << iter << " CPU " << cpu << "\n";
        
        if ( precondCounter == N_TEST )
        {
            // if the differential of times is significant, use the fastest method
            // but otherwise, select the simplest method:
            precondMethod = 0;
            if ( precondCPU[0] > precondCPU[1] + 10 )
                precondMethod = 1;
            if ( precondCPU[precondMethod] > precondCPU[2] + 10 )
                precondMethod = 2;
            
            if ( prop->verbose )
            {
                std::clog << " precond 0 time " << precondCPU[0] << "\n";
                std::clog << "         1 time " << precondCPU[1] << "\n";
                std::clog << "         2 time " << precondCPU[2] << "\n";
                std::clog << " ----> " << precondMethod << std::endl;
            }
        }
        else
        {
            //alternate betwen methods { 0, 1, 2 }
            precondMethod = ( 1 + precondMethod ) % 3;
        }
    }
    else if ( precondCounter > PERIOD )
    {
        precondCPU[0] = 0;
        precondCPU[1] = 0;
        precondCPU[2] = 0;
        precondCPU[3] = 0;
        precondCounter = 0;
    }
}


void Simul::apply() const
{
    sMeca.apply();
}


void Simul::computeForces() const
{
    // we could use here an accessory Meca mec;
    try {
        prop->complete(*this);
        setInteractions(sMeca);
        sMeca.computeForces();
    }
    catch ( Exception & e )
    {
        std::clog << "cytosim could not compute the forces:\n";
        std::clog << "   " << e.what() << std::endl;
    }
}


void Simul::dump() const
{
    std::string cwd = FilePath::get_cwd();
    const char path[] = "dump";
    FilePath::make_dir(path);
    FilePath::change_dir(path);
    sMeca.dump();
    FilePath::change_dir(cwd);
    std::clog << "cytosim dumped a system of size " << sMeca.dimension() << "in " << path << std::endl;
}


//==============================================================================
//                              SOLVE-X 1D
//==============================================================================


void Simul::solveX()
{
    //-----initialize-----
    
    sMeca1D.clear();
    
    for(Fiber * fib = fibers.first(); fib; fib=fib->next())
        sMeca1D.add(fib);

    sMeca1D.prepare(prop->time_step, prop->kT);
    
    //-----set matrix-----

    for ( Couple * co = couples.firstAA(); co ; co=co->next() )
    {
        Hand const* h1 = co->hand1();
        Hand const* h2 = co->hand2();
        
        const index_t i1 = h1->fiber()->matIndex();
        const index_t i2 = h2->fiber()->matIndex();
        assert_true( i1 != i2 );
        
        sMeca1D.addLink(i1, i2, co->stiffness(), h2->pos().XX - h1->pos().XX);
    }
    
    for ( Single * gh = singles.firstA(); gh ; gh=gh->next() )
    {
        Hand const* h = gh->hand();
        const index_t ii = h->fiber()->matIndex();
        
        sMeca1D.addClamp(ii, gh->prop->stiffness, gh->position().XX - h->pos().XX);
    }
    
    //-----resolution-----

    real noise = sMeca1D.setRightHandSide(prop->kT);
    
    sMeca1D.solve(prop->tolerance * noise);
    sMeca1D.apply();
}
