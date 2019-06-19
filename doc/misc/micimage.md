# Pseudo Microscope Images
 
 The accessory program `micimage` can be used to display a system as if it was imaged in a fluoresence light microscope.
 
 Producing the fake micrographs involves 3 steps:
 - running the simulation with `sim`
 - generating files containing coordinates of 'fluorophores' with `reportF`
 - render these coordinate files as grayscale PNG images with `micimage`
 .

 Example:
 
	 sim pombe.cym
	 reportF fiber:speckles interval=0.1
	 
 `reportF` will generate files 'report????.txt' containing the positions of the speckles
 The value of `interval` dertermines the average distance between speckles, and thus the
 linear density of labelling along the fiber.
 The speckles are distributed randomly and uniformly along the filaments.
 It is possible then to convert these coordinate files to images:
 
 	micimage report00*txt dim=256,128
 
 Options to `micimage` can be used to adjust the result, by changing the magnification, etc.
 
	micimage help
 
 You can then @ref DisplayMovie "assemble the images into a movie".
 