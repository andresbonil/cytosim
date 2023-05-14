# Cytosim - Microtubule Breaking Functionality 

![Cytosim](doc/data/cytosim.png)

This branch of Cytosim contains microtubule breaking functionality based on cultured cell in-vitro experiments from the Verhey Lab at Michigan Medicine's Cell and Developmental Biology department. Live imaging data has revealed that clustered kinesin can cause microtubule crosslinking, sliding, and breakage. Through this plugin, we hope to obtain data that influences future experimental design, along with utilizing experimental data and literature to fine-tune representations of the microtubule breaking phenomena. 

# Fiber.cc Added Functionality 

We have created an algorithm which utilizes a user-selected threshold in order to determine the probability of certain segments in microtubules to break within a simulation. During each frame of a simulation, when breaking is enabled, tensions at segments are evaluated, and if that tension value passes the probability check, a severing event is queued and executed. These changes are located in ```/src/sim/fiber.cc```.

# Changes in Fiber Objects

Fiber props have been given two new parameters, ```breaking``` and ```breaking_threshold```. ```breaking``` is a boolean value that specifies whether breaking is enabled for said fiber object (0 for false, 1 for true). ```breaking_threshold``` is an int that is the force threshold (in piconewtons) needed for a microtubule to break. These parameters can be specified for any fiber object in each simulation's config.cym file. The addition of these parameters can be found in ```fiber_prop.h```/```fiber_prop.cc```.

An example of specifying afiber object with breaking enabled and a threshold of 20 piconewtons in config.cym would be as follows:

set fiber microtubule
	{
	    rigidity = 10
	    segmentation = 0.1
	    breaking = 1
		breaking_threshold = 20
	}

# Changes in Report Executable

Through testing this plugin, we have found it useful to output the total number of fibers in each frame of the simulation to assess the quantity of breaking events. Thus, we have added the function ```reportFiberNum()``` to ```simul_report.cc``` and ```simul.h```. Through executing ```./report fiber:num``` and parsing the output, we can receieve informative output relating to microtubule quantity over time. 


# Contributors
This plugin was developed by Andres Bonilla and Qi (Archie) Geng with assistance from Dr. Kristen Verhey and Dr. Francois Nedelec 

# Contact
Email: abonil@umich.edu

# Cytosim Documentation/Troubleshooting 
Please refer to the original release of Cytosim by Dr. Francois Nedelec for documentation regarding the application and installation: https://gitlab.com/f-nedelec/cytosim
