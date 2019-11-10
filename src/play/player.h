// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef PLAYER_H
#define PLAYER_H

///this turns on display code in simul.h
#define  DISPLAY

#include "simul.h"
#include "parser.h"
#include "display.h"
#include "sim_thread.h"
#include "display_prop.h"
#include "play_prop.h"

class FiberDisp;
class View;

/// Base of the Graphical User Interface for Cytosim
class Player
{
public:

    /// the display parameters
    DisplayProp  DP;
    
    /// the parameters for play
    PlayProp     PP;
    
    /// container for object's display properties
    PropertyList dproperties;
    
    /// SimThread to control the live simulation
    SimThread    thread;
    
    /// a flag for live simulation
    bool         goLive;

    /// an alias to `thread` Simul
    Simul    &   simul;
    
    /// the Display object
    Display  *   mDisplay;
    
    //---------------------------------COMMANDS---------------------------------
    
    /// constructor
    Player();

    ///
    ~Player();
    
    /// initialize display
    void initialize();
  
    /// cleanup
    void clear();

    /// return all Fiber FiberDisp
    PropertyList allFiberDisp();
   
    /// return all Fiber FiberDisp
    PropertyList allVisibleFiberDisp();
 
    /// return all Hand PointDisp
    PropertyList allHandDisp();
    
    /// return all Hand PointDisp for which 'visible==true'
    PropertyList allVisibleHandDisp();

    /// return all Sphere/Solid/Bead PointDisp
    PropertyList allSphereDisp();
    
    /// return all Space PointDisp
    PropertyList allSpaceDisp();
    
    /// return a FiberDisp
    FiberDisp * firstFiberDisp();

    //---------------------------------COMMANDS---------------------------------
    
    /// reset view, without changing the current frame
    void rewind();
   
    /// start animation
    bool startPlayback();
    
    /// accelerate animation if possible
    void accelerate();
    
    /// start reverse animation
    bool startBackward();
    
    /// start live simulation
    void extendLive();

    /// stop animation
    void stop();

    /// start or stop animation
    void startstop();
    
    /// reset the sim-state and timer
    void restart();

    /// load previous frame
    void previousFrame();
    
    /// go to the next frame, returns 1 if EOF is reached
    void nextFrame();
    
    /// write global display parameters
    void writePlayParameters(std::ostream& out, bool prune) const;

    /// write Object display parameters
    void writeDisplayParameters(std::ostream& out, bool prune) const;
    
    //-----------------------------DISPLAY--------------------------------------
  
    /// initialize display with given style
    void setStyle(int);

    /// build message that appears on top
    std::string buildReport(std::string) const;
    
    /// build message that appears on bottom
    std::string buildLabel() const;
    
    /// build central message
    std::string buildMemo(int) const;

    
    /// set View::focus and quat to match the center of gravity of the Fibers
    void autoTrack(FiberSet const&, View&);
    
    /// adjust the viewing area
    void autoScale(SpaceSet const&, View&);

    /// adjust the model view and load frame if asked
    void prepareDisplay(View&, int mag);
    
    /// read parameters contained in string
    void readDisplayString(View&, std::string const&);
    
    /// display cytosim state and message
    void displayCytosim();
    
    /// display function calls displayCytosim
    void displayScene(View&, int mag);
    
    /// export current viewport to a graphic file
    int  saveView(const char * root, unsigned indx, int verbose = 1) const;
    
    /// save high-resolution image using composite method
    int  saveViewMagnified(int mag, const char* filename, const char* format, int downsample=1);
    
    /// save a high-resolution image using composite method
    int  saveViewMagnified(int mag, const char* root, unsigned indx, int downsample=1);

};


#endif
