# Stochastic events in Cytosim
 
Most stochastic events are simulated as such, using pseudo random numbers generated on the fly.
 
If a stochastic event occurs at a constant rate, its time of occurence in generated as:

	time = -log(random()) / rate

where `random()` returns a random number uniformly distributed in ]0,1] (zero is excluded). This is the 'inverse transform sampling' method to generate exponential variates. In this case, the variable `time` is exponentially distributed with expectancy `1/rate`.


# Kramers reaction rate theory

The detachment of a molecular link follows the theory described by Hendrik Kramers in:

>  Brownian motion in a field of force and the diffusion model of chemical reactions  
>  H.A. Kramers - Physica VII, no 4, pp284-304 - 1940

Essentially, the detachment rate `off_rate` varies with the force exerted on the link:
 
	off_rate = unbinding_rate * exp( force / unbinding_force )
 
Where `force` is the norm of the force vector calculated by `cytosim`,
while `unbinding_rate` and `unbinding_force` are constants associated with the bound state.
In cytosim, these parameters are specified in the properties of `Hands`.

# Time-dependent rates

Using Kramers rate theory implies that the rate of the event is varying in time.  
The Gillespie approach needs to be modified, and we follow the procedure described in:
 
>  A Dynamical Monte Carlo Algorithm for Master Equations with Time-Dependent Transition Rates  
>  A. Prados et al. Journal of Statistical Physics, Vol. 89, Nos. 3/4, 1997  
>  http://dx.doi.org/10.1007/BF02765541

 
In short, a normalized time `esp` is first generated,
again using a random number uniformly distributed in [0,1] provided by `random()`.
At each time step, `esp` is reduced as a function of the value of the rate during the interval.
The associated event is performed if `esp` becomes negative, which is when the associated time
crosses the 'present' into the past.
  
	time = 0;
	esp = -log(random());
	
	while ( time < max_time )
	{
		time += time_step;
		esp -= time_step * rate(time);
		while ( esp < 0 )
		{
		    do_event();
		    esp += -log(random())
		}
	}
 
In this example, the event can be performed multiple times in the same time_step, but this would not be done for detachment and other events that can only occur once. 
