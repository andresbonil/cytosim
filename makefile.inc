# Cytosim was created by Francois Nedelec. Copyright 2007-2018 EMBL.


#---------------- MACHINE = {mac, linux, cluster, cygwin, auto}
MACHINE := auto

#---------------- COMPILER = {gcc, icc, clang}
COMPILER := gcc

#---------------- MODE = {F = Fast; D = Debug; P = Profiling; C = Coverage}
MODE := D

#---------------- INTEL Math Kernel Library
# If the Intel Math Kernel library is a free library providing BLAS+LAPACK,
# If it is installed on your system, you can set:
#     HAS_MKL := 1 for dynamic linking
#     HAS_MKL := 2 for static linking
# otherwise, set HAS_MKL := 0;

HAS_MKL := 0

#---------------- PNG image export support
# `libpng` needs to be installed to save PNG images:
#     Mac OSX:
#        www.macports.org :  port install libpng 
#        Homebrew  (http://brew.sh)
#             brew install libpng
#             brew tap homebrew/dupes
#             brew install zlib
#     CentOS:  yum install libpng-devel
#     Ubuntu:  apt-get or synaptic libpng-dev
# HAS_PNG can be:
#     0 : no PNG support
#     1 : macports installation
#     2 : Homebrew installation

HAS_PNG := 0

#-------------------------------------------------------------------------------
#---------------------------  Platform Detection  ------------------------------
#-------------------------------------------------------------------------------

ifeq ($(MACHINE), auto)

    ifeq ($(shell uname), Darwin)
       MACHINE := mac
    else ifeq ($(shell uname -n), login.cluster.embl.de)
       MACHINE := cluster
    else ifeq ($(shell uname), Linux)
       MACHINE := linux
    else ifeq ($(shell uname -o), Cygwin)
       MACHINE := cygwin
    else
       $(error Unknown platform)
    endif

endif

#-------------------------------------------------------------------------------
#------------------------  Compiler Specifications  ----------------------------
#-------------------------------------------------------------------------------

# common options for all C++ objects:
CXXFLG := -std=gnu++11 -fno-rtti -Wredeclared-class-member -Wno-deprecated-declarations
CXXFLG := -std=gnu++11 -fno-rtti -fopenmp #-g -fno-inline
CXXFLG := -std=gnu++11 -fno-rtti

# option concerning memory alignment:
ALIGN  := -malign-double -falign-loops -falign-jumps -falign-functions
ALIGN  :=

# options concerning warnings:
WARN   := -Wundef -Wall -Wno-unknown-pragmas# -Wno-unused

ifeq ($(MACHINE),mac)
    LIBTOOL := /usr/bin/libtool -static -o
    CXXFLG := $(CXXFLG) -fno-pic
else
    LIBTOOL := libtool --mode=link gcc -g -O -o
    LIBTOOL := ar -rcs
endif


ifeq ($(COMPILER),gcc)

# Intel SIMD: -march=core2 (SSE3); -march=sandybridge (AVX); -march=haswell (AVX2)
# You can use g++-8 to support Intel OMP directives
    CXX      := g++ -Wfatal-errors -Wno-deprecated-declarations
    FAST     :=-O3 -fno-trapping-math -ffast-math -funroll-loops #-floop-optimize2 -ftree-vectorize -mavx2 -mfma
    FAST     :=-O3 -fno-trapping-math -ffast-math -funroll-loops -march=native -ffp-contract=fast

    DEBUG    :=-g3 -ggdb
    COVERAGE :=-fprofile-arcs -ftest-coverage
    
else ifeq ($(COMPILER),icc)

    LIBTOOL  := ar -rcs

    #CXX     := icpc -fp-trap=common -vec-report0
    #CXX     := icpc -bind_at_load -fvisibility=hidden -wd1224,161
    CXX      := icpc -fp-trap=none -vec-report0 -wd1224,161,1478
    CXX      := icpc -fargument-noalias -fp-trap=none -wd1224,161,1478
    # Intel advanced instruction sets:
    # '-xHost' to optimize for host machine
    # '-xAVX' for AVX
    # '-march=corei7-avx' for AVX
    # '-march=core-avx2' for Intel core i7 (ca. 2015)
    FAST     := -O3 -march=core-avx2
    FAST     := -O3 -march=core-avx2 -Wl,-no_pie -funroll-loops -fno-math-errno -fp-model fast=1 -fp-model no-except
    WARN     := -qopt-report=2 -qopt-report-phase=vec# -opt-report-file=stderr
    WARN     :=

    DEBUG    := -march=core-avx2 -g -wd279,383,810,869,981,1418,1419,1572,2259
    COVERAGE :=

else ifeq ($(COMPILER),clang)

# for the standard library, use: -stdlib=libc++
# the old library: -stdlib=libstdc++ (deprecated)

    CXX      := c++ -Wfatal-errors
    FAST     := -O3 -ffast-math -funroll-loops -mavx2 -mfma
    FAST     := -O3 -ffast-math -funroll-loops
    DEBUG    := -g3 -ggdb
    COVERAGE := -fprofile-arcs -ftest-coverage

endif

#-------------------------------------------------------------------------------
#------------------------------- Mac OSX ---------------------------------------
#-------------------------------------------------------------------------------

ifeq ($(MACHINE),mac)

    MKLDIR := /opt/intel/mkl/lib

    FlagsC := $(COVERAGE)
    FlagsD := $(DEBUG) $(WARN)
    FlagsP := -O2 $(ALIGN) -fno-inline
    FlagsF := $(FAST) $(ALIGN)

    LINK   := -framework Accelerate

    LINKGL := -framework GLUT -framework OpenGL -framework AGL
    # Using brew's freeGLUT:
    #LINKGL = -L/usr/local/Cellar/freeglut/3.0.0 -lglut -framework OpenGL
    
    ifeq ($(HAS_PNG), 1)

        # macports libraries:
        LIB_PNG := /opt/local/lib/libpng.a /opt/local/lib/libz.a
        INC_PNG := -I/opt/local/include

    endif
	
    ifeq ($(HAS_PNG), 2)

        # HomeBrew libraries:
        LIB_PNG := /usr/local/lib/libpng.a /usr/local/opt/zlib/lib/libz.a
        INC_PNG := -I/usr/local/include

    endif

endif

#-------------------------------------------------------------------------------
#--------------------------------- Cluster -------------------------------------
#-------------------------------------------------------------------------------

ifeq ($(MACHINE),cluster)

    # MKLROOT should be defined by 'module load imkl'
    MKLDIR := $(MKLROOT)/lib/intel64

    # Beware of AVX2 and other advanced features on heterogeneous clusters:
    FAST   :=-O3 -fno-tree-vectorize -ffast-math -funroll-loops -mavx -ffp-contract=fast

    FlagsD := -O0 $(DEBUG) $(WARN)
    FlagsP := $(FAST) $(ALIGN) -O0 -fno-inline -pg
    FlagsC := -O0 -fprofile-arcs -ftest-coverage
    FlagsF := $(FAST) $(ALIGN) -finline-functions
    
    ### self-compiled BLAS/LAPACK fortran library:
    # this works on the new cluster (64 bit, gcc 4.1 or higher)
    #LINK  := -static -L/g/nedelec/opt/netlib/linux_em64t -llapack -lblas -lgfortran
      
    ### 64bit Linux distribution, 2011 (distribution on LSF EMBL cluster)
    # ln -sf /usr/lib64/liblapack.so.3 /usr/lib64/liblapack.so
    # ln -sf /usr/lib64/libblas.so.3 /usr/lib64/libblas.so
    #LINK  := /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3 -L/usr/lib64 -lpthread
    
    ### SLURM cluster (EMBL, 2017)
    LINK   := -llapack -lopenblas -lgfortran -lpthread
    ### SLURM cluster with self-compiled LAPACK/BLAS (EMBL, 2018)
    LINK   := /g/nedelec/lapack/liblapack.a /g/nedelec/lapack/librefblas.a -lgfortran -lpthread

    # linkage options for programs using graphics:
    #LINKGL: =-L/usr/lib -L/usr/lib64 -lglut -lGL -lGLU -lXt -lX11

endif

#-------------------------------------------------------------------------------
#---------------------------------- Linux --------------------------------------
#-------------------------------------------------------------------------------

ifeq ($(MACHINE),linux)

    ### Linux options
    FlagsD := -O0 $(DEBUG) $(WARN)
    FlagsP := -O0 -fno-inline -pg $(ALIGN) $(FAST)
    FlagsC := -O0 -fprofile-arcs -ftest-coverage
    FlagsF := -O3 -finline-functions $(ALIGN) $(FAST)

    # specify the lib path:
    USRLIB := /usr/lib

    ### For 64bit Centos 5.8 March 2012:
    # Centos  blas: yum install blas-devel
    # Centos  lapack: yum install lapack-devel
    ### 64bit Linux distribution Ubuntu 12.04 LTS, 2012
    # Ubuntu   blas : apt-get install libblas.dev
    # Ubuntu lapack : apt-get install liblapack.dev

    #### For dynamic linking:
    LINK := -L$(USRLIB) -llapack -lblas -lpthread
    #### For static linking specify the objects directly:
    LINK := -L$(USRLIB) $(USRLIB)/liblapack.so $(USRLIB)/libblas.so -lpthread

    ### Install libraries needed for graphical executables:
    # Centos freglut: yum install freeglut-devel
    # Centos libXi:   yum install libXi-devel
    # Centos libXmu:  yum install libXmu-devel

    # Ubuntu freglut: apt-get install freeglut3-dev
    # Ubuntu libXi:   apt-get install libXi-dev
    # Ubuntu libXmu:  apt-get install libXmu-dev

    ### linkage options for programs using graphics:
    #   check for correct library libname.so otherwise give full path
    LINKGL := -L$(USRLIB) -lglut -lGL -lGLU -lXt -lX11 -lGLEW

    ifneq ($(HAS_PNG), 0)

        LIB_PNG := $(USRLIB)/libpng.a $(USRLIB)/libz.a
        ING_PNG :=

    endif

endif


#-------------------------------------------------------------------------------
#------------------------------- Cygwin ----------------------------------------
#-------------------------------------------------------------------------------

ifeq ($(MACHINE),cygwin)

    ### Cygwin October 2016
    FlagsD := -O0 $(DEBUG) $(WARN)
    FlagsP := -O0 -fno-inline -pg $(ALIGN) $(FAST)
    FlagsC := -O0 -fprofile-arcs -ftest-coverage
    FlagsF := -O3 -finline-functions $(ALIGN) $(FAST)

    ### Needs
    # blas:   yum install blas-devel
    # lapack: yum install lapack-devel
    LINK := -L/usr/lib -L/bin -llapack -lblas -lpthread -lgfortran

    ### linkage options for programs using graphics:
    ### Needs
    ### freglut: yum install freeglut-devel
    ### libXi:   yum install libXi-devel
    ### libXmu:  yum install libXmu-devel
    ###  check for correct library yourlibname.so otherwise give full path
    LINKGL := -L/usr/lib -lglut -lGL -lGLU -lXt -lX11 -lGLEW

    ifneq ($(HAS_PNG), 0)

    LIB_PNG := /lib/libpng.dll.a /lib/libz.a
    ING_PNG :=

    endif

    BINEXT :=.exe

endif


#-------------------------------------------------------------------------------
#------------------------------- Native Windows --------------------------------
#-------------------------------------------------------------------------------
# THIS IS NOT MAINTAINED: it worked in 2008 but has not been tried since

ifeq ($(MACHINE),windows)

    # Use a windows-native blas/lapack compilation:
    # copying "libblas.a" and "liblapack.a" in "/usr/lib"

    FlagsD := -g  $(WARN)
    FlagsP := -pg -O3 -march=i686 -ffast-math $(ALIGN)
    FlagsF := -O3 -finline-functions -march=i686 -ffast-math $(ALIGN)
  
    LINK   := -llapack -lblas
    LINKGL := -lglut32 -lglu32 -lopengl32
    BINEXT :=.exe

endif



#-------------------------------------------------------------------------------
#------------------ linking with Intel Math Kernel Library ---------------------
#-------------------------------------------------------------------------------

# The Intel MKL provides LAPACK and BLAS implementations optimized for Intel chips


ifeq ($(HAS_MKL),1)

    # sequential dynamic linking:
    MKLLIB := -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
    
    # threaded dynamic linking:
    #MKLLIB := -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

    # modify the linking command:
    LINK := -L$(MKLDIR) $(MKLLIB) -lpthread

endif


ifeq ($(HAS_MKL),2)

    # This are options for the compiler:
    MKLOPT := -m64 -I${MKLROOT}/include

    # sequential static linking (newer GCC):
    MKLLIB := -Wl,--start-group $(MKLDIR)/libmkl_intel_lp64.a $(MKLDIR)/libmkl_sequential.a $(MKLDIR)/libmkl_core.a -Wl,--end-group

    # sequential static linking (older GCC):
    MKLLIB := $(MKLDIR)/libmkl_intel_lp64.a $(MKLDIR)/libmkl_sequential.a $(MKLDIR)/libmkl_core.a

    # threaded static linking:
    #MKLLIB := $(MKLDIR)/libmkl_intel_lp64.a $(MKLDIR)/libmkl_thread.a $(MKLDIR)/libmkl_core.a /opt/intel/lib/libiomp5.a

    # modify the linking command:
    LINK := $(MKLLIB) -lpthread -lm -ldl

endif

