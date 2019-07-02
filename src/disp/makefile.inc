# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#
# File src/disp/makefile.inc

OBJ_DISPLAY=gle.o gle_color.o gle_color_list.o view.o view_prop.o glapp.o


#-----------------------LIB & DEF for PNG support-------------------------------

IMAGE_LIB=
IMAGE_DEF=

ifneq ($(HAS_PNG), 0)

    IMAGE_DEF+=-DHAS_PNG $(INC_PNG)
    IMAGE_LIB+=$(LIB_PNG)

endif


#----------------------------targets--------------------------------------------


libcytodisp.a: $(OBJ_DISPLAY) offscreen.o saveimage.o grid_display.o
	$(MAKELIB)
	$(DONE)


saveimage.o: saveimage.cc saveimage.h
	$(COMPILE) $(IMAGE_DEF) -c $< -o build/$@


offscreen.o: offscreen.cc offscreen.h offscreen_fbo.cc offscreen_glx.cc
	$(COMPILE) -c $< -o build/$@


$(OBJ_DISPLAY) grid_display.o: %.o: %.cc %.h | build
	$(COMPILE) $(addprefix -Isrc/, math base) -c $< -o build/$@