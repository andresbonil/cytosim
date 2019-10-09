# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


TESTS:=test test_gillespie test_solve test_random test_math test_glos test_quaternion\
       test_code test_matrix test_thread test_blas test_pipe

TESTS_GL:=test_opengl test_vbo test_glut test_glapp test_platonic\
          test_rasterizer test_space test_grid test_sphere


.PHONY: tests
tests: $(TESTS) $(TESTS_GL)

vpath %.cc src/test


#--------------------macros----------------------------------------------------

TEST_MAKE = $(COMPILE) $(addprefix -Isrc/, math base sim disp play) $(OBJECTS) $(LINK) -o bin/$@

TEST_GL = -I src/disp $(GRAFIX) $(IMAGES)


#----------------------------Targets--------------------------------------------

test: test.cc random.o SFMT.o tictoc.o backtrace.o
	$(TEST_MAKE)
	$(DONE)
vpath test bin

test_blas: test_blas.cc random.o SFMT.o backtrace.o
	$(TEST_MAKE)
	$(DONE)
vpath test_blas bin

test_code: test_code.cc random.o SFMT.o tictoc.o 
	$(COMPILE) $(INC_BM) $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)
vpath test_code bin

test_matrix: test_matrix.cc matsparsesymblk.o matsparsesym.o matsparsesym1.o matsparsesym2.o matrix.o random.o SFMT.o tictoc.o backtrace.o
	$(TEST_MAKE)
	$(DONE)
vpath test_matrix bin

test_glos: test_glos.cc glossary.o filepath.o tokenizer.o stream_func.o exceptions.o backtrace.o ansi_colors.o
	$(TEST_MAKE)
	$(DONE)
vpath test_glos bin

test_quaternion: test_quaternion.cc cytomath.a SFMT.o cytobase.a 
	$(TEST_MAKE)
	$(DONE)
vpath test_quaternion bin

test_random: test_random.cc random.o SFMT.o exceptions.o filewrapper.o messages.o smath.o backtrace.o tictoc.o
	$(TEST_MAKE)
	$(DONE)
vpath test_random bin

test_asm: test_asm.cc
	$(COMPILE) -restrict -S -c -m64 $^ -o build/test_asm.s
	$(DONE)
vpath test_asm.s build

test_simd: test_simd.cc random.o SFMT.o tictoc.o
	$(COMPILE) -Isrc/base -Isrc/math -Isrc/SFMT $(OBJECTS) -o bin/$@
	$(DONE)
vpath test_simd bin

test_omp: test_omp.cc
	$(COMPILE) -fopenmp $^ -o bin/test_omp
	$(DONE)
vpath test_omp bin

test_sizeof: test_sizeof.cc
	$(TEST_MAKE)
	$(DONE)
vpath test_sizeof bin

test_pipe: test_pipe.cc gle_color.cc gle_color_list.cc
	$(TEST_MAKE)
	$(DONE)
vpath test_pipe bin

test_cxx: test_cxx.cc
	$(TEST_MAKE)
	$(DONE)
vpath test_cxx bin

test_math: test_math.cc
	$(TEST_MAKE)
	$(DONE)
vpath test_math bin

test_thread: test_thread.cc
	$(TEST_MAKE)
	$(DONE)
vpath test_thread bin

test_string: test_string.cc
	$(TEST_MAKE)
	$(DONE)
vpath test_string bin

test_gillespie: test_gillespie.cc random.o SFMT.o backtrace.o
	$(TEST_MAKE)
	$(DONE)
vpath test_gillespie bin

#----------------------------Targets with Graphics------------------------------

OBJGLAP := $(OBJ_DISPLAY) offscreen.o

test_opengl: test_opengl.cc
	$(CXX) $^ $(GRAFIX) -o bin/$@
	$(DONE)
vpath test_opengl bin

test_vbo: test_vbo.cc
	$(CXX) $^ $(GRAFIX) -o bin/$@
	$(DONE)
vpath test_vbo bin

test_glut: test_glut.cc
	$(TEST_MAKE) $(TEST_GL)
	$(DONE)
vpath test_glut bin

test_glapp: test_glapp.cc $(OBJGLAP) cytobase.a cytomath.a
	$(TEST_MAKE) $(TEST_GL)
	$(DONE)
vpath test_glapp bin

test_sphere: test_sphere.cc $(OBJGLAP) cytomath.a cytobase.a
	$(TEST_MAKE) $(TEST_GL)
	$(DONE)
vpath test_sphere bin

test_platonic: test_platonic.cc platonic.o $(OBJGLAP) cytomath.a cytobase.a
	$(TEST_MAKE) $(TEST_GL)
	$(DONE)
vpath test_platonic bin

test_rasterizer: test_rasterizer.cc rasterizerG.o cytomath.a $(OBJGLAP) cytobase.a
	$(TEST_MAKE) $(TEST_GL)
	$(DONE)
vpath test_rasterizer bin

test_space: test_space.cc cytosim.a cytomath.a $(OBJGLAP) cytobase.a
	$(TEST_MAKE) $(TEST_GL)
	$(DONE)
vpath test_space bin

test_grid: test_grid.cc tictoc.o cytomath.a cytobase.a $(OBJGLAP) grid_display.o
	$(TEST_MAKE) $(TEST_GL)
	$(DONE)
vpath test_grid bin

test_solve: test_solve.cc cytomath.a cytobase.a
	$(TEST_MAKE) $(TEST_GL)
	$(DONE)
vpath test_solve bin

