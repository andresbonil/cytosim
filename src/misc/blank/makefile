# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
# makefile  Francois Nedelec 2007.

#-------------------include compiler-specifications-----------------------------

include makefile.inc

#-------------------check compiler and compiler-flags---------------------------

ifeq ($CXX,)
   $(error No compiler defined for $$(MACHINE)=$(MACHINE))
endif

ifeq ($(Flags$(MODE)),)
   $(error No compiler-options defined for $$(MACHINE)=$(MACHINE) in mode $(MODE))
endif

COMPILE := $(CXX) $(Flags$(MODE))

#----------------------------objects--------------------------------------------

INCDIR=import base disp

vpath %.h  $(INCDIR)
vpath %.cc $(INCDIR)
vpath %.c  $(INCDIR)
vpath %.o  build

OBJ=platonic.o random.o exceptions.o vector3.o\
    glossary.o filepath.o parameter.o iomessages.o filewrapper.o\
    gle.o gle_color.o glapp.o object.o blank.o

#----------------------------targets--------------------------------------------

blank: $(OBJ) SFMT.o
	$(COMPILE) $^ $(GRAFIX) -o $@
vpath blank .

.PHONY: clean
clean:
	rm -f blank build/*

.PHONY: tar
tar: clean
	tar -cf blank.tar.bz2 --dereference * --bzip2 


#----------------------------rules----------------------------------------------

INCLUDE=$(addprefix -I, $(INCDIR))

$(OBJ): %.o: %.cc
	$(COMPILE) $(INCLUDE) -c $< -o build/$@

SFMT.o: SFMT.c SFMT-std.c SFMT-sse2.c SFMT-alti.c
	$(CXX) -x c -DNDEBUG -DSFMT_MEXP=19937 -Iimport -c $< -o build/$@

#----------------------------makedep--------------------------------------------

MAKEDEP=gcc -MM
SRC=$(patsubst %.o,%.cc, $(OBJ))  #$(addsuffix .h, $(BASE))
FILES=$(findfile $(SRC), $(INCDIR))
FILES=$(findfile random.h, $(INCDIR))

build/makefile.dep: $(SRC)
	$(MAKEDEP) $(INCLUDE) $(filter %.cc, $^) >  $@

-include build/makefile.dep
