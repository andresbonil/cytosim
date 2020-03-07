# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#
# File src/disp/makefile.inc

OBJ_DISP := fiber_disp.o line_disp.o point_disp.o grid_display.o\
            gle.o gle_color.o gle_color_list.o view.o view_prop.o glapp.o


#-----------------------LIB & DEF for PNG support-------------------------------

IMAGES=
IMAGES_DEF=

ifneq ($(HAS_PNG), 0)

    IMAGES_DEF+=-DHAS_PNG $(INC_PNG)
    IMAGES+=$(LIB_PNG)

endif


#----------------------------targets--------------------------------------------


cytodisp.a: $(OBJ_DISP) saveimage.o offscreen.o
	$(MAKELIB)
	$(DONE)


saveimage.o: saveimage.cc saveimage.h
	$(COMPILE) $(IMAGES_DEF) -c $< -o build/$@


offscreen.o: offscreen.cc offscreen.h offscreen_fbo.cc offscreen_glx.cc
	$(COMPILE) -c $< -o build/$@


$(OBJ_DISP): %.o: %.cc %.h | build
	$(COMPILE) -Isrc/math -Isrc/base -Isrc/sim -c $< -o build/$@
