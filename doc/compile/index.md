# How to compile Cytosim

The core of Cytosim is written in C++ and it is necessary to recompile the programs after each modification of the source code. Some of the accessory tools use [Python](https://www.python.org).
 
Compilation requires a C++ compiler: e.g. [gcc](http://gcc.gnu.org/), [Clang](http://clang.llvm.org) or the [Intel compiler](http://en.wikipedia.org/wiki/Intel_C%2B%2B_Compiler), together with a few libraries.
Compilation is started from a terminal, with a program called [make](http://www.gnu.org/software/make/).

### Dimensionality

The dimensionality is backed into the executable during compilation. 
To change it, [follow these instructions](dimensionality.md).

### Mathematical libraries
 
Running cytosim with `sim` or `play` requires these mathematical libraries:

* [BLAS](http://netlib.org/blas)
* [LAPACK](http://netlib.org/lapack)

These libraries offer standard interface to linear algebra functions used in Cytosim.
There is a [public reference implementations](http://netlib.org). 
Compiling the reference code is possible with a [FORTRAN](http://en.wikipedia.org/wiki/Fortran) compiler.
Precompiled libraries are available for most platforms, within:

- [Intel Math Kernel Library](http://software.intel.com/en-us/articles/intel-mkl/)
- [Apple's vecLib](http://developer.apple.com/hardwaredrivers/ve/vector_libraries.html)
- [OpenBLAS](https://www.openblas.net)
- Also available for many Linux distributions.

Apple's veclib is preinstalled on Mac OSX, and [Intel's MKL](https://en.wikipedia.org/wiki/Math_Kernel_Library) is available free of charge. MKL has better performance on Intel processors.

### Graphical libraries
 
Cytosim's ***play*** relies on [OpenGL](http://www.opengl.org/) for display 
and uses [POSIX threads](http://en.wikipedia.org/wiki/POSIX_Threads) for multithreading.
 
Interactivity is provided by [GLUT](http://www.opengl.org/resources/libraries/glut/),
which can be replaced by [freeGLUT](http://freeglut.sourceforge.net/).
 
GLUT and OpenGL are included in Mac OSX:

	-framework GLUT -framework OpenGL 

Cytosim program `play` will be able to export images only if it is linked with the PNG graphical libraries, but everything else works fine. See the [compile options](options.md).

# Getting Ready 

### Mac-OSX

On Mac OSX, install [Xcode](https://developer.apple.com/technologies/tools/), which is available on the Mac App Store. After installing Xcode, install the Xcode 'Command-Line Tools', which is an optional package providing 'make'. All necessary libraries are already included in MacOSX.

We provide the Xcode project file for cytosim, which is a convenient way to access the code.

Optionally, Cytosim can use the mouse wheel to zoom in and out, if you use FreeGLUT, instead of Apple's GLUT, or [Renaud Blanch's GLUT patch](http://iihm.imag.fr/blanch/howtos/MacOSXGLUTMouseWheel.html).


### Linux

On Linux, you need to install the GNU compiler collection, BLAS and LAPACK. This is sufficient to compile `sim`.
Recent Linux distributions provide precompiled BLAS/LAPACK as optional installation (check your distribution).

To make `play` install the OpenGL developer libraries and FreeGLUT or OpenGLUT. 


### Windows

Native compilation on Windows is a Herculean task. We thus recommend installing Cygwin, which emulates a Linux-like environment. A fresh installation of Cygwin is recommended to ensure that all required reference libraries are installed.  
You will need a compiler, the X window system, BLAS/LAPACK and GLUT.  
Please, [refer to the dedicated page](cygwin.md).


# Compilation

After installing a compiler and `gnu's make`, 
you are ready to compile from a terminal, with the following commands in the root directory of cytosim:

	make sim
	make play

The command `make` without arguments will build `sim` and `play`.  
You can then check the resulting executables, that should be located in subdirectory `bin`:

	bin/sim info
	bin/sim
	bin/play live

# Troubleshooting 1

Compilation is specified in `makefile` and `makefile.inc`, and these files may need to be adjusted manually.
Check  `makefile.inc` first and verify the [compile options](options.md).

It attepts to automatically detect the platforms:

	#---------------- MACHINE = {mac, linux, cygwin, auto}
	MACHINE := auto

Manually set `MACHINE` to `mac`, `linux` or `cygwin` depending on your platform,
and check the parameters set lower in the `makefile.inc`, for this platform, for example:

	ifeq ($(MACHINE),linux)
		...
	endif

To poinpoint the problem, try to build the objects in the order of how many depencies they have:

#### Check your compiler for C++11 support and compilation switches

    make test_cxx

#### Check your BLAS/LAPACK directly:

	make test_blas
	
#### Check if it can read a configuration file:

	make test_glos

#### At this stage, you can compile `sim`:

	make sim

# Troubleshooting 2: `play`

If you are having trouble compiling `play`, check its requirements independently:

#### Check for thread support:

	make test_thread
	
#### Check for OpenGL support:

	make test_opengl
	
#### Check for GLUT support:

	make test_glut

#### Check for GLAP support:

	make test_glap

#### At this stage, you can compile `play`:

	make play

### Finally: contact us!

Please, write to `feedbackATcytosimDOTorg`.

Please, describe what fails and what you have tried.
Attach your 'makefile.inc' and tell us the platform on which you compiled.
