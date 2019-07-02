# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


TOOLS:=frametool sieve reader report reportF

.PHONY: tools
tools: $(TOOLS)

vpath %.cc src/tools

#--------------------macros----------------------------------------------------

TOOL_OBJ := cytosim.a cytospace.a cytomath.a matsparsesymblk.o cytobase.a

TOOL_INC  = $(addprefix -Isrc/, math base sim sim/spaces sim/fibers sim/singles sim/organizers disp play)

TOOL_MAKE = $(COMPILE) $(TOOL_INC) $(OBJECTS) $(LINK) -o bin/$@


#-------------------targets----------------------------------------------------
 
 
frametool: frametool.cc
	$(COMPILE) $^ -o bin/$@
	$(DONE)
vpath frametool bin


sieve: sieve.cc frame_reader.o $(TOOL_OBJ)
	$(TOOL_MAKE)
	$(DONE)
vpath sieve bin


reader: reader.cc frame_reader.o $(TOOL_OBJ)
	$(TOOL_MAKE)
	$(DONE)
vpath reader bin


report: report.cc frame_reader.o $(TOOL_OBJ) | bin
	$(COMPILE) $(TOOL_INC) $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)
vpath report bin


bin1/report: report.cc frame_reader.o cytosimD1.a cytospaceD1.a cytomathD1.a cytobase.a
	$(COMPILE) -DDIM=1 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin1/report
	$(DONE)
vpath bin1/report bin1


bin2/report: report.cc frame_reader.o cytosimD2.a cytospaceD2.a cytomathD2.a cytobase.a
	$(COMPILE) -DDIM=2 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin2/report
	$(DONE)
vpath bin2/report bin2


bin3/report: report.cc frame_reader.o cytosimD3.a cytospaceD3.a cytomathD3.a cytobase.a
	$(COMPILE) -DDIM=3 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin3/report
	$(DONE)
vpath bin3/report bin3


reportF: reportF.cc frame_reader.o $(TOOL_OBJ)
	$(TOOL_MAKE)
	$(DONE)
vpath reportF bin


cymart: cymart.cc frame_reader.o $(TOOL_OBJ)
	$(TOOL_MAKE)
	$(DONE)
vpath cymart bin


#----------------------------makedep--------------------------------------------

dep/part7.dep:
	$(MAKEDEP) $(wildcard src/tools/*.cc) > $@
