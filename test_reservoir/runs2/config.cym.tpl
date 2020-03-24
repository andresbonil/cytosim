

set simul system
{
    time_step = 0.005
    viscosity = 0.1
    precondition = 0
    display = ( style=2; )
}

set space cell
{
    shape = circle
}

set fiber microtubule
{
    rigidity = 20
    segmentation = 0.5
    lattice = 1, 0.008
}

set hand hand1
{
    binding_rate = 9
    binding_range = 0.1
    unbinding_rate = 0.2

    activity = walk
    step_size = 0.008
	footprint = 0
    unloaded_speed = 0
	stall_force=1

    hold_growing_ends = 0.99

    display = ( color=green; size=7; width=7; )
}

set hand hand2
{
    binding_rate = 9
    binding_range = 0.1
    unbinding_rate = 0.2

    activity = walk
    step_size = 0.008
	footprint = 0
    unloaded_speed = 0
	stall_force=1
    hold_growing_ends = 0.99

    display = ( color=light_blue; size=7; width=7; )
}

set single single1
{
    hand = hand1
    diffusion = fast
    reservoir_add = 10
}

set single single2
{
    hand = hand2
    diffusion = fast
    reservoir_add = 10
}

new cell
{
    radius = 10
}

new 1 microtubule
{
    length = 10
    orientation = 1
    position = 0
}

new 10 single1
new 10 single2

run 4000 system
{
    nb_frames = 100
    solve = 0
}

