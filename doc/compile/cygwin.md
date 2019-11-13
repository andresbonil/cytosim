## Instructions to run Cytosim on Cygwin

Cytosim runs on Windows, using [Cygwin](https://www.cygwin.com), but this is not supported directly by us.

### Installation

Run the [Cygwin installer](https://cygwin.com/install.html) with default settings until you get to the package selection screen. To compile and run Cytosim, several packages must be explicitely selected:

- gcc				(GNU compiler suite including C++)
- make				(GNU make)
- libBLAS			(Basic Linear Algebra Subprograms)
- libLAPACK		(Linear Algebra PACKage)
- OpenGL			(Open Graphics Library)
- GLEW				(OpenGL Extension Wrangler Library)
- GLUT				(OpenGL Utility Toolkit)
- X11    	    	(X Window System)
- python3     	(python interpreter for cygwin)

Search for the keywords provided and select one matching item mentionning
"devel" or "lib". To select items for installation click the weird icon to the left of "Skip" or "Default" until it reads "Install" instead.
After all packages are selected press 'Next' to install Cygwin.  

**Tip**: In the view drop-down menu select "Category" which will reorganize packages to show them by category. You may then chose to install all the developer tools for Cygwin, which is overkill and will take ~20-30 minutes, but may be easier than finding out the tools listed above.

### Compilation

Compile from within the Cygwin terminal, which will use the toolchain (`gcc` and `make`) provided by cygwin. 
The MACHINE should be automatically selected (in case of trouble, set `MACHINE:=cygwin` in `makefile.inc`). Disable offscreen rendering by setting `HAS_PNG:=0`.The procedure is the same as on other platforms (enter `make`). If you experience trouble, please let us know.

### The X Window System (X11)

Finally, all graphical tools included in cytosim (in particular 'play') will use the [X Window System](https://en.wikipedia.org/wiki/X_Window_System) also known as X11 to open a window, and this will work only if a X11 server is running on your local computer. Two X11 implementations can be installed with Cygwin:

[Cygwin/X](https://en.wikipedia.org/wiki/Cygwin/X)

[Xming](https://en.wikipedia.org/wiki/Xming)

We got things to work with the free 2007 version of Xming. Xming can be installed while Cygwin installation is running. 
Once Xming is installed run XLaunch to start the X11 window server, and direct your application to this server by setting the `DISPLAY` environment variable.

This will instruct 'play' to send X11 queries to the default X11 server on the local machine. This command should be run once every time a new cygwin terminal is opened, otherwise 'play' will fail with a message "failed to open display". 

**Tip**: To debug X11, it may be useful to run standard program such as `xclok`. This way you can test if the X11 server is operational.

### Starting

Cytosim will only work within Cygwin, and not directly under Windows. This means that it cannot be started from the Windows command line, and you can close the Windows/MS-DOS command prompt. 

Cytosim should be started from the Cygwin terminal window:

    cd cytosim
    export DISPLAY=:0
    bin/play live   

Author: Daniel Cortes and Francois Nedelec, last updated 7 Nov 2019
