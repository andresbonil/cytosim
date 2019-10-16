// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "frame_reader.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"


// Use the second definition to get some verbose reports:
#define VLOG(ARG) ((void) 0)
//#define VLOG(ARG) std::clog << ARG;

//------------------------------------------------------------------------------

FrameReader::FrameReader() : inputter(DIM)
{
    frameIndex = -1;
}

void FrameReader::clear()
{
    inputter.rewind();
    clearPositions();
}

void FrameReader::openFile(std::string const& file)
{
    clearPositions();
    
    int error = inputter.open(file.c_str(), "rb");
    
    if ( error )
    {
        //file was not found, we try 'gunzip'
        std::string tmp = file + ".gz";
        FILE* fp = fopen(tmp.c_str(), "r");
        if ( fp )
        {
            fclose(fp);
            tmp = "gunzip " + tmp;
            std::clog << tmp << std::endl;
            
            if ( 0 == system(tmp.c_str()) )
                inputter.open(file.c_str(), "rb");
        }
    }
    
    if ( !inputter.file() )
        throw InvalidIO("file `"+file+"' not found");

    if ( inputter.error() )
        throw InvalidIO("file `"+file+"' is invalid");
 
    inputter.vectorSize(DIM);
    //std::clog << "FrameReader: has openned " << obj_file << std::endl;
}


int FrameReader::badFile()
{
    if ( !inputter.file() )
        return 8;
    
    if ( inputter.eof() )
        inputter.clear();
    
    if ( ! inputter.good() )
        return 7;
    
    return 0;
}


void FrameReader::checkFile()
{
    if ( !inputter.file() )
        throw InvalidIO("No open file");
    
    if ( inputter.eof() )
        inputter.clear();
    
    if ( ! inputter.good() )
        throw InvalidIO("File has errors");
}

//------------------------------------------------------------------------------
#pragma mark -

void FrameReader::clearPositions()
{
    VLOG("FrameReader: clear\n");
    
    frameIndex = -1;
    framePos.clear();
    framePos.reserve(1024);
}


void FrameReader::savePos(int frm, const fpos_t& pos, int s)
{
    if ( frm < 0 )
        return;
    
    size_t inx = frm;
    
    if ( inx >= framePos.capacity() )
    {
        constexpr size_t chunk = 1024;
        size_t sz = ( inx + chunk - 1 ) & ~( chunk -1 );
        framePos.reserve(sz);
    }
    
    if ( inx >= framePos.size() )
    {
        size_t i = framePos.size();
        framePos.resize(inx+1);
        while ( i <= inx )
            framePos[i++].status = 0;
    }
    
    if ( framePos[inx].status < s )
    {
        framePos[inx].status = s;
        framePos[inx].position = pos;
    
        //VLOG("FrameReader: position of frame " << frm << " is " << pos << '\n');
        VLOG("FrameReader: learned position of frame " << frm << '\n');
    }
}


/**
 This uses the current knowledge to move to a position
 in the file where we should find frame `frm`.
*/
int FrameReader::seekPos(int frm)
{
    if ( inputter.eof() )
        inputter.clear();
    
    if ( frm < 1 || framePos.empty() )
    {
        VLOG("FrameReader: seekPos rewind\n");
        inputter.rewind();
        return 0;
    }
    
    int inx = std::min(frm, lastPossibleFrame());

    while ( inx > 0  &&  framePos[inx].status == 0 )
        --inx;
    
    //check if we know already were the frame starts:
    if ( 0 < inx )
    {
        VLOG("FrameReader: using known position of frame " << inx << '\n');
        inputter.set_pos(framePos[inx].position);
        return inx;
    }
    else {
        VLOG("FrameReader: rewind\n");
        inputter.rewind();
        return 0;
    }
}


int FrameReader::lastKnownFrame() const
{
    int res = lastPossibleFrame();
    while ( 0 < res  &&  framePos[res].status < 2 )
        --res;
    return res;
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 scan file forward from current position to find the next Cytosim frame
 @return 0 if no frame was found
*/
int FrameReader::seekFrame(const int frm)
{        
    VLOG("FrameReader: seekFrame("<< frm <<")\n");
    
    int inx = seekPos(frm);
    
    if ( inx == frm )
        return 0;
    
    while ( ! inputter.eof() )
    {
        fpos_t pos;
        bool has_pos = false;
        std::string line;

        do {
            has_pos = !inputter.get_pos(pos);
            line = inputter.get_line();
            
            if ( inputter.eof() )
                return 1;
            
#ifdef BACKWARD_COMPATIBILITY // 2012
            if ( 0 == line.compare(0, 7, "#frame ") )
                break;
#endif
            
        } while ( line.compare(0, 9, "#Cytosim ") );
        
        //std::clog << "******\n";
        VLOG("FrameReader: " << line << '\n');

        if ( ! inputter.eof() )
        {
            if ( has_pos ) savePos(inx, pos, 2);
            if ( inx == frm )
            {
                if ( has_pos ) inputter.set_pos(pos);
                return 0;
            }
            ++inx;
        }
    }
    
    VLOG("FrameReader: seekFrame("<< frm <<") reached EOF\n");
    return 1;
}

//------------------------------------------------------------------------------
/** 
 returns 0 for success, an error code, or throws an exception
 */
int FrameReader::loadFrame(Simul& sim, int frm, const bool reload)
{
    if ( badFile() )
        return 7;

    VLOG("FrameReader: loadFrame("<<frm<<", " << reload <<")\n");
    
    // a negative index is counted from the last frame
    if ( frm < 0 )
    {
        int res = loadLastFrame(sim);
        if ( frm == -1 ) return res;
        VLOG("FrameReader: counting down from frame " << lastKnownFrame() << '\n');
        frm = std::max(0, frm + 1 + lastKnownFrame());
    }
    
    // what we are looking for might already be in the buffer:
    if ( frm == frameIndex && ! reload )
        return 0;
    
    // it might be the next one in the buffer:
    if ( frm == 1+frameIndex )
        return loadNextFrame(sim);

    // otherwise, try to find the start tag from there:
    if ( 0 != seekFrame(frm) )
        return 1;
    
    // store the position in the file:
    fpos_t pos;
    bool has_pos = !inputter.get_pos(pos);
    
    VLOG("FrameReader: reading frame " << frm << " from " << pos << '\n');
    //VLOG("FrameReader: reading frame " << frm << '\n');
    
    // ask cytosim to read the file:
    if ( 0 == sim.reloadObjects(inputter) )
    {
        VLOG("FrameReader: loadFrame("<< frm <<") successful\n");
        frameIndex = frm;
        if ( has_pos ) savePos(frameIndex, pos, 3);

        // the next frame should start at the current position:
        has_pos = !inputter.get_pos(pos);
        if ( has_pos ) savePos(frameIndex+1, pos, 1);
        return 0;
    }
    else
    {
        VLOG("FrameReader: loadFrame("<< frm <<") EOF at frame " << frm << '\n');
        return 1;
    }
}


/** 
 returns 0 for success, an error code, or throws an exception
 */
int FrameReader::loadNextFrame(Simul& sim)
{
    if ( badFile() )
        return 7;
    
    fpos_t pos;
    bool has_pos = !inputter.get_pos(pos);

    if ( 0 == sim.reloadObjects(inputter) )
    {
        ++frameIndex;
        
        // the position we used was good, to read this frame
        if ( has_pos ) savePos(frameIndex, pos, 3);

        VLOG("FrameReader: loadNextFrame() after frame " << currentFrame() << '\n');
        
        // the next frame should start from the current position:
        if ( !inputter.get_pos(pos) ) savePos(frameIndex+1, pos, 1);
        return 0;
    } 
    else
    {
        VLOG("FrameReader: loadNextFrame() EOF after frame " << currentFrame() << '\n');
        return 1;
    }
}

/**
 returns 0 for success, an error code, or throws an exception
 */
int FrameReader::loadLastFrame(Simul& sim)
{
    if ( badFile() )
        return 7;
    
    /// seek last known position:
    int frm = lastKnownFrame();
    if ( frm > 0 )
        inputter.set_pos(framePos[frm].position);
    else
        inputter.rewind();
    
    int res = 1;
    while ( 0 == loadNextFrame(sim) )
        res = 0;
    
    return res;
}
