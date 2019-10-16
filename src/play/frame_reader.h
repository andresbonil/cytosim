// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "iowrapper.h"
#include <vector>

class Simul;


/// Helper class to access a particular frame in a trajectory file
/** 
 FrameReader is used to find a particular frame (eg. frame 10) in a trajectory file.
 FrameReader will handle basic IO failures and will remember the starting points 
 of all frames that were sucessfully found, using this information to speed up 
 future access to these and other frames in the same file. For example, if the 
 position of frame 50 is known and frame 60 is requested, it will start searching
 the file from the end of frame 50.

 FrameReader makes minimal assuptions on what constitutes a 'frame':
 - It looks for a string that identifies the beggining of a frame (Cytosim).
 - It calls Simul::reloadObjects() to actually read the content of the frame.
 .
*/
class FrameReader
{    
private:
    
    /// accessory class to store a position in a file
    class file_pos 
    {
    public:
        int    status;   ///< indicates if `position` is valid
        fpos_t position; ///< starting position in the file
        file_pos() { status=0; }
    };
    
    /// type for list of positions
    typedef  std::vector<file_pos> PosList;
    
    /// the stream from which input is made
    Inputter inputter;
    
    /// starting position for each frame
    PosList  framePos;
    
    /// index of frame stored currently
    int      frameIndex;

    /// remember position `pos` as the place where frame `frm` should start
    void     savePos(int frm, const fpos_t& pos, int status);
   
    /// go to a position where a frame close to `frm` is known to start
    int      seekPos(int frm);
    
    /// check file validity
    void     checkFile();
    
    /// return 0 if file is good for input
    int      badFile();
    
    /// return last
    int      lastPossibleFrame() const { return framePos.size()-1; }
    
public:
    
    /// constructor, after which openFile() should be called
    FrameReader();
    
    /// open file for input
    void     openFile(std::string const& file);
    
    /// clear the buffer
    void     clearPositions();
    
    /// last frame seen in the file
    int      lastKnownFrame() const;
    
    /// return state of file object
    bool     hasFile() { return inputter.file(); }
    
    /// true when end of file is reached
    bool     eof() const  { return inputter.eof();  }
    
    /// rewind file
    void     rewind() { inputter.rewind(); frameIndex=-1; }

    /// true if everything looks correct for input
    bool     good() const { return inputter.good(); }

    /// rewind file and clear position buffer
    void     clear();
    
    /// return index of current frame 
    long     currentFrame() const { return frameIndex; }

    /// true if current frame is valid
    bool     hasFrame() const { return frameIndex >= 0; }
    
    /// find the starting point of frame `frm` by brute force !
    int      seekFrame(int frm);
    
    /// load specified frame into Simul (frm=0: first; frm=-1: last)
    int      loadFrame(Simul&, int frm, bool reload = false);
    
    /// read the next frame in the file, return 1 for EOF
    int      loadNextFrame(Simul&);
    
    /// read the last frame in the file, return 1 if no frame was found
    int      loadLastFrame(Simul&);

};


