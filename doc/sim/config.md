# Cytosim's Configuration File

The virtual system and the sequence of actions to be performed is specified in a configuration file. 
This config file must be a plain (ASCII) text file, as produced by many text editor (but not word processors).
By default, it is called **`config.cym`**.
 
Cytosim understand a [small set of commands](commands.md) and [predefined objects](objects.md) with their [associated parameters](parameters.md).

One needs to call at least:

 - `set` to define a new object category, and set its parameter values,
 - `new` to create objects,
 - `run` to perform simulation steps. 

About the config file:

- Parameters should be specified in [units of `(s, um, pN)`](units.md)
- Parameter names are predefined in the simulations for each [object class](objects.md)
- Parameters are specified with a ` = ` sign, and can contain multiple values
- Two parameters can be specified on the same line if they are separated by a semi-column: `;`
- Curly brackets `{  }` are used to group the parameters together in logical units
- Parentheses `(  )` can be used for subgroups, such as `display`
- Strings can be encapsulated by double quotes `" "` or parentheses  

# Editing

Cytosim's configuration files `*.cym` are plain text ASCII file, and should be edited using a [PLAIN TEXT editor](https://en.wikipedia.org/wiki/Text_editor).
Do not use word processors such as Microsoft Word. 

Recommended editors:

- [TextMate](https://macromates.com)
- [Atom](https://atom.io)
- [Gleany](https://www.geany.org)
- [Sublime text](https://www.sublimetext.com)

Manuel Lera-Ramirez made a [syntax highlighting configuration](../misc/Cytosim.tmbundle.zip) for TextMate.

# Example
 
Many examples can be found in the directory [***cym/***](../../cym).
	 

	% Self organization of Microtubules driven by bivalent Motors
	% Adapted from Nedelec et al. Nature, 1998
	
	set simul system
	{
		time_step = 0.01
		viscosity = 0.05
		display = ( style=2; )
	}
	
	set space cell
	{
		shape = circle
	}
	
	new cell
	{   
		radius = 10
	}

	set fiber microtubule
	{
		rigidity = 20
		segmentation = 0.5
		confine = inside, 100
		display = ( line_width=1; color=white; )
	}
	
	set hand kinesin
	{
		binding = 10, 0.01   % rate, range
		unbinding = 0.1, 3   % rate, force
		
		activity = move
		unloaded_speed = 0.8
		stall_force = 5
	
		bind_also_end = 1
		hold_growing_end = 1
	
		display = ( size=7; width=7; )
	}
	
	set couple complex
	{
		hand1 = kinesin
		hand2 = kinesin
		stiffness = 100
		diffusion = 10
	}
	
	new 100 microtubule
	{
		length = 9
	}
	
	new 2000 complex
	
	run 5000 system
	{
		nb_frames = 50
	}
	
