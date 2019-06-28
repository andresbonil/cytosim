# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


OBJ_SPACE=space.o space_prop.o space_square.o space_sphere.o space_disc.o\
          space_dice.o space_torus.o space_polygon.o space_polygonZ.o \
          space_ellipse.o space_cylinder.o space_ring.o space_cylinderZ.o space_capsule.o \
          space_strip.o space_periodic.o space_banana.o space_cylinderP.o\

OBJ_SPACED1=$(OBJ_SPACE:.o=D1.o)
OBJ_SPACED2=$(OBJ_SPACE:.o=D2.o)
OBJ_SPACED3=$(OBJ_SPACE:.o=D3.o)

OBJ_SPACEG=$(OBJ_SPACE:.o=G.o)
OBJ_SPACED1G=$(OBJ_SPACE:.o=D1G.o)
OBJ_SPACED2G=$(OBJ_SPACE:.o=D2G.o)
OBJ_SPACED3G=$(OBJ_SPACE:.o=D3G.o)

#----------------------------targets--------------------------------------------

cytospace.a: $(OBJ_SPACE)
	$(MAKELIB)
	$(DONE)


cytospaceD1.a: $(OBJ_SPACED1)
	$(MAKELIB)
	$(DONE)


cytospaceD2.a: $(OBJ_SPACED2)
	$(MAKELIB)
	$(DONE)


cytospaceD3.a: $(OBJ_SPACED3)
	$(MAKELIB)
	$(DONE)


#library with display
cytospaceG.a: $(OBJ_SPACEG)
	$(MAKELIB)
	$(DONE)


#1D library with display
cytospaceD1G.a: $(OBJ_SPACED1G)
	$(MAKELIB)
	$(DONE)


#2D library with display
cytospaceD2G.a: $(OBJ_SPACED2G)
	$(MAKELIB)
	$(DONE)


#3D library with display
cytospaceD3G.a: $(OBJ_SPACED3G)
	$(MAKELIB)
	$(DONE)


#----------------------------rules----------------------------------------------
INC_BMD=$(addprefix -Isrc/, base math sim sim/spaces disp)

$(OBJ_SPACE): %.o: %.cc %.h | build
	$(COMPILE) $(INC_BMD) -c $< -o build/$@

$(OBJ_SPACED1): %D1.o: %.cc %.h | build
	$(COMPILE) -DDIM=1 $(INC_BMD) -c $< -o build/$@

$(OBJ_SPACED2): %D2.o: %.cc %.h | build
	$(COMPILE) -DDIM=2 $(INC_BMD) -c $< -o build/$@

$(OBJ_SPACED3): %D3.o: %.cc %.h | build
	$(COMPILE) -DDIM=3 $(INC_BMD) -c $< -o build/$@


$(OBJ_SPACEG): %G.o: %.cc %.h | build
	$(COMPILE) -DDISPLAY $(INC_BMD) -c $< -o build/$@

$(OBJ_SPACED1G): %D1G.o: %.cc %.h | build
	$(COMPILE) -DDIM=1 -DDISPLAY $(INC_BMD) -c $< -o build/$@

$(OBJ_SPACED2G): %D2G.o: %.cc %.h | build
	$(COMPILE) -DDIM=2 -DDISPLAY $(INC_BMD) -c $< -o build/$@

$(OBJ_SPACED3G): %D3G.o: %.cc %.h | build
	$(COMPILE) -DDIM=3 -DDISPLAY $(INC_BMD) -c $< -o build/$@

