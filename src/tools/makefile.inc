# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


TOOLS:=frametool sieve reader report reportF report_python

.PHONY: tools
tools: $(TOOLS)

vpath %.cc src/tools

#--------------------macros----------------------------------------------------

TOOL_OBJ := cytosim.a cytomath.a matsparsesymblk.o cytobase.a SFMT.o

TOOL_INC  = $(addprefix -Isrc/, math base sim sim/spaces sim/fibers sim/singles sim/organizers disp play)

TOOL_MAKE = $(COMPILE) $(TOOL_INC) $(OBJECTS) $(LINK) -o bin/$@

#PYTHONI = -fPIC -shared  -I/home/dmitrief/conda3/include/python3.7m -I/home/dmitrief/.local/lib/python3.7/site-packages/pybind11/include

TEST := $(shell python3 -m pybind11 --includes)
PYTHONI := -fPIC -shared $(TEST) -Wl,-undefined,dynamic_lookup
SUFFIX := $(shell python3-config --extension-suffix)
#-------------------targets----------------------------------------------------

frametool: frametool.cc | bin
	$(COMPILE) $^ -o bin/$@
	$(DONE)
vpath frametool bin


sieve: sieve.cc frame_reader.o $(TOOL_OBJ) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath sieve bin


reader: reader.cc frame_reader.o $(TOOL_OBJ) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath reader bin


report: report.cc frame_reader.o $(TOOL_OBJ) | bin
	$(COMPILE) $(TOOL_INC) $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)
vpath report bin

report_python: report_python.cc frame_reader.o sim_thread.o $(TOOL_OBJ) | bin
	$(COMPILE) $(PYTHONI) $(TOOL_INC) $(OBJECTS) $(LINK) -fPIC $(INCS) -o bin/cytosim$(SUFFIX)
	$(DONE)
vpath report_python bin


bin1/report: report.cc frame_reader.o cytosimD1.a cytomathD1.a cytobase.a SFMT.o | bin1
	$(COMPILE) -DDIM=1 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin1/report
	$(DONE)
vpath bin1/report bin1


bin2/report: report.cc frame_reader.o cytosimD2.a cytomathD2.a cytobase.a SFMT.o | bin2
	$(COMPILE) -DDIM=2 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin2/report
	$(DONE)
vpath bin2/report bin2


bin3/report: report.cc frame_reader.o cytosimD3.a cytomathD3.a cytobase.a SFMT.o | bin3
	$(COMPILE) -DDIM=3 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin3/report
	$(DONE)
vpath bin3/report bin3


reportF: reportF.cc frame_reader.o $(TOOL_OBJ) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath reportF bin


cymart: cymart.cc frame_reader.o $(TOOL_OBJ) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath cymart bin


#----------------------------makedep--------------------------------------------

dep/part7.dep:
	$(MAKEDEP) $(wildcard src/tools/*.cc) > $@

