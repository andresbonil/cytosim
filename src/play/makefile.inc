# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

OBJ_PLAY := fiber_disp.o line_disp.o point_disp.o display_prop.o play_prop.o\
            sim_thread.o display.o display1.o display2.o display3.o frame_reader.o

FILES_PLAYER := player.cc player.h player_disp.cc

OBJ_PLAY1D := $(OBJ_PLAY:.o=D1.o)
OBJ_PLAY2D := $(OBJ_PLAY:.o=D2.o)
OBJ_PLAY3D := $(OBJ_PLAY:.o=D3.o)

PLAY_INC = $(addprefix -Isrc/, math base sim disp play sim/organizers sim/singles)
PLAY_DEP = $(addprefix -Isrc/, math base sim disp play) $(LINK) $(LINKGL) $(IMAGE_LIB)

################################ objects with DIM defined in source code:

player.o: $(FILES_PLAYER) | build
	$(COMPILE) $(PLAY_INC) -c $< -o build/$@

$(OBJ_PLAY): %.o: %.cc %.h | build
	$(COMPILE) $(PLAY_INC) -c $< -o build/$@

play: play.cc player.o $(OBJ_PLAY) cytosim.a cytospaceG.a libcytodisp.a cytomath.a cytobase.a | bin
	$(COMPILE) $(INFO) $(OBJECTS) $(PLAY_DEP) -o bin/play
	$(DONE)
vpath play bin

################################ objects with DIM = 1:

playerD1.o: $(FILES_PLAYER) | build
	$(COMPILE) -DDIM=1 $(PLAY_INC) -c $< -o build/$@

$(OBJ_PLAY1D): %D1.o: %.cc %.h | build
	$(COMPILE) -DDIM=1 $(PLAY_INC) -c $< -o build/$@

bin1/play: play.cc playerD1.o $(OBJ_PLAY1D) cytosimD1.a cytospaceD1G.a libcytodisp.a cytomathD1.a cytobase.a
	$(COMPILE) -DDIM=1 $(INFO) $(OBJECTS) $(PLAY_DEP) -o bin1/play
	$(DONE)
vpath bin1/play bin1

################################ objects with DIM = 2:

playerD2.o: $(FILES_PLAYER) | build
	$(COMPILE) -DDIM=2 $(PLAY_INC) -c $< -o build/$@

$(OBJ_PLAY2D): %D2.o: %.cc %.h | build
	$(COMPILE) -DDIM=2 $(PLAY_INC) -c $< -o build/$@

bin2/play: play.cc playerD2.o $(OBJ_PLAY2D) cytosimD2.a cytospaceD2G.a libcytodisp.a cytomathD2.a cytobase.a
	$(COMPILE) -DDIM=2 $(INFO) $(OBJECTS) $(PLAY_DEP) -o bin2/play
	$(DONE)
vpath bin2/play bin2

################################ objects with DIM = 3:

playerD3.o: $(FILES_PLAYER) | build
	$(COMPILE) -DDIM=3 $(PLAY_INC) -c $< -o build/$@

$(OBJ_PLAY3D): %D3.o: %.cc %.h | build
	$(COMPILE) -DDIM=3 $(PLAY_INC) -c $< -o build/$@

bin3/play: play.cc playerD3.o $(OBJ_PLAY3D) cytosimD3.a cytospaceD3G.a libcytodisp.a cytomathD3.a cytobase.a
	$(COMPILE) -DDIM=3 $(INFO) $(OBJECTS) $(PLAY_DEP) -o bin3/play
	$(DONE)
vpath bin3/play bin3
