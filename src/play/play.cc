// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "simul_prop.h"
#include "glossary.h"
#include "messages.h"
#include "offscreen.h"
#include "saveimage.h"
#include "filepath.h"
#include "splash.h"
#include "ansi_colors.h"

#include "opengl.h"
#include "player.h"
#include "view.h"
#include "gle.h"

Player player;

SimThread& thread = player.thread;
Simul&      simul = thread.sim();
PlayProp&      PP = player.PP;
DisplayProp&   DP = player.DP;

/// enable to create a player for command-line-only offscreen rendering
//#define HEADLESS_PLAYER

#ifdef HEADLESS_PLAYER
void helpKeys(std::ostream& os) {}
#else
#  include "glut.h"
#  include "glapp.h"
#  include "fiber_prop.h"
#  include "fiber_disp.h"
#  include "point_disp.h"
using glApp::flashText;
void buildMenus();
#  include "play_keys.cc"
#  include "play_menus.cc"
#  include "play_mouse.cc"
#endif

void goodbye()
{
    //printf("Goodbye...\n");
    player.clear();
}

//------------------------------------------------------------------------------
#pragma mark - Display

/**
 display is done only if data can be accessed by current thread
 */
void displayLive(View& view, int mag)
{
    if ( 0 == thread.trylock() )
    {
        // read and execute commands from incoming pipe:
        thread.readInput(32);
        //thread.debug("display locked");
        if ( simul.prop->display_fresh )
        {
            player.readDisplayString(view, simul.prop->display);
            simul.prop->display_fresh = false;
        }
        //thread.debug("display");
        player.prepareDisplay(view, mag);
        player.displayCytosim();
        thread.unlock();
        glFinish();
    }
    else
    {
        //thread.debug("display: trylock failed");
        //glutPostRedisplay();
    }
}


/**
 This is a bare-bone version used for off-screen rendering.
 */
void displayOffscreen(View & view, int mag)
{
    //std::clog << "displayOffscreen " << glApp::views.size() << '\n';
    player.displayScene(view, mag);
    glFinish();
}


/// copy data from multisample to normal buffer
void blitBuffers(GLuint normal, GLuint multi, GLint W, GLint H)
{
    //std::clog << "blitting multisample buffer\n";
    glBindFramebuffer(GL_READ_FRAMEBUFFER, multi);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, normal);
    glBlitFramebuffer(0, 0, W, H, 0, 0, W, H, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, multi);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, normal);
}

//------------------------------------------------------------------------------
#pragma mark - main


void help(std::ostream& os)
{
    os << "play [OPTIONS] [PATH] [FILE]\n"
          "     live                     enter live simulation mode directly\n"
          "     PATH                     change working directory as specified\n"
          "     FILE.cym                 specify input configuration file\n"
          "     FILE.cmo                 specify trajectory file\n"
          "     FILE.cyp                 specify display configuration file\n"
          "     PARAMETER=value          set parameter value (example size=512)\n"
          "     image frame=INT          render specified frame offscreen\n"
          "     image frame=INT,INT,...  render several frames offscreen\n"
          "     image magnify=INT        render frames at higher resolution\n"
          "     movie                    render all frames off screen\n"
          "     movie=on                 render all frames on screen\n"
          "     movie period=INT         render one frame every INT frames\n"
          " (there should be no whitespace around the equal sign)\n";
}


/// different modes:
enum { ONSCREEN, OFFSCREEN_IMAGE, OFFSCREEN_MOVIE, ONSCREEN_MOVIE };


int main(int argc, char* argv[])
{
    int  mode = ONSCREEN;
    int  magnify = 1;
    Glossary arg;
    
    Cytosim::all_silent();
    std::atexit(goodbye);

    if ( arg.read_strings(argc-1, argv+1) )
        return EXIT_FAILURE;
    
    // check for major options:
    
    if ( arg.use_key("help") )
    {
        splash(std::cout);
        help(std::cout);
        return EXIT_SUCCESS;
    }

    if ( arg.use_key("info") || arg.use_key("--version") )
    {
        splash(std::cout);
        print_version(std::cout);
        if ( SaveImage::supported("png") )
            std::cout << "    PNG enabled\n";
        return EXIT_SUCCESS;
    }
    
    if ( arg.use_key("live") || arg.has_key(".cym") )
        player.goLive = true;
    
    if ( arg.use_key("image") )
        mode = OFFSCREEN_IMAGE;

    if ( arg.use_key("poster") )
    {
        mode = OFFSCREEN_IMAGE;
        magnify = 3;
    }
    
    if ( arg.value_is("movie", 0, "on") )
        mode = ONSCREEN_MOVIE;
    else if ( arg.use_key("movie") )
        mode = OFFSCREEN_MOVIE;
    
    // get image over-sampling:
    arg.set(magnify, "magnify") || arg.set(magnify, "magnification");

    // change working directory if specified:
    if ( arg.has_key("directory") )
    {
        FilePath::change_dir(arg.value("directory", 0));
        //std::clog << "Cytosim working directory is " << FilePath::get_cwd() << '\n';
    }

    // can specify a frame index to be loaded:
    size_t frm = 0;
    bool has_frame = false;
    
    try
    {
        has_frame = arg.set(frm, "frame");
        simul.initialize(arg);
    }
    catch( Exception & e )
    {
        print_magenta(std::cerr, "Error: "+e.brief());
        std::cerr << e.info() << '\n';
        return EXIT_FAILURE;
    }
    
    // default configuration file for play:
    std::string setup = player.goLive ? simul.prop->config_file : simul.prop->property_file;
    
#ifdef HEADLESS_PLAYER
    View view("*");
    view.setDisplayFunc(displayOffscreen);
#else
    glApp::setDimensionality(DIM);
    if ( arg.use_key("fullscreen") )
        glApp::setFullScreen(1);
    View& view = glApp::views[0];
    view.setDisplayFunc(displayLive);
#endif

    try
    {
        // check for play's configuration file specified on the command line:
        bool has_setup = arg.set(setup, ".cyp");
        
        // extract first specification of "simul:display" string from the setup file
        if ( FilePath::is_file(setup) )
            Parser(simul, 0, 0, 0, 0, 0).readConfig(setup);
        else if ( has_setup )
            std::cerr << " warning: could not read `" << setup << "'\n";
        
        // read settings from the setup file, but do not overwrite the command-line options:
        arg.read(simul.prop->display, 1);
        simul.prop->display_fresh = false;
        
        if ( !arg.empty() )
        {
            view.read(arg);
            DP.read(arg);
            PP.read(arg);
        }
    }
    catch( Exception & e )
    {
        print_magenta(std::cerr, "Error: "+e.brief());
        std::cerr << e.info() << '\n';
        return EXIT_FAILURE;
    }
    
    //---------Open trajectory file and read state

    if ( ! player.goLive || has_frame )
    {
        try
        {
            std::string file = simul.prop->property_file;
            
            if ( !FilePath::is_file(file) )
                throw InvalidIO("could not find `"+file+"'\n");
            
            Parser(simul, 1, 1, 0, 0, 0).readConfig(file);
            
            // read 'setup' file again allowing to overwrite 'display' values
            if ( file != setup )
                Parser(simul, 0, 0, 0, 0, 0).readConfig(setup);
            
            thread.openFile(simul.prop->trajectory_file);
            
            // load requested frame from trajectory file:
            if ( thread.loadFrame(frm) )
            {
                // load last frame in file:
                if ( thread.loadLastFrame() )
                    std::cerr << "Warning: could only load frame " << thread.currentFrame() << ' ';
            }
            frm = thread.currentFrame();
        }
        catch( Exception & e )
        {
            arg.warnings(std::cerr);
            print_magenta(std::cerr, "Error: "+e.brief());
            std::cerr << e.info() << '\n';
            return EXIT_FAILURE;
        }
    }
    
#ifndef __APPLE__
    // it is necessary under Linux/Windows to initialize GLUT to display fonts
    glutInit(&argc, argv);
#endif
    
    //-------- off-screen (non interactive) rendering -------
    
    if ( mode == OFFSCREEN_IMAGE || mode == OFFSCREEN_MOVIE )
    {
        const int W = view.width() * magnify;
        const int H = view.height() * magnify;
        
        //std::cerr << W << "x" << H << << " downsample " << PP.downsample << '\n';
        
        if ( !OffScreen::openContext() )
        {
            std::cerr << "Failed to create off-screen context\n";
            return EXIT_FAILURE;
        }
        GLuint fbo = OffScreen::createBuffer(W, H, 0);
        if ( !fbo )
        {
            std::cerr << "Failed to create off-screen pixels\n";
            return EXIT_FAILURE;
        }
        GLuint multi = 0;
        if ( view.multisample > 1 )
        {
            multi = OffScreen::createBuffer(W, H, view.multisample);
        }
        
        gle::initialize();
        player.setStyle(DP.style);
        view.initGL();

        if ( mode == OFFSCREEN_IMAGE )
        {
            unsigned inx = 0;
            // it is possible to specify multiple frame indices:
            do {
                thread.loadFrame(frm);
                // only save requested frames:
                if ( thread.currentFrame() == frm )
                {
                    displayOffscreen(view, magnify);
                    if ( multi )
                        blitBuffers(fbo, multi, W, H);
                    if ( magnify > 1 )
                        player.saveView("poster", frm, 2);
                    else
                        player.saveView("image", frm, 2);
                }
            } while ( arg.set(frm, "frame", ++inx) );
        }
        else if ( mode == OFFSCREEN_MOVIE )
        {
            // save every PP.period
            unsigned s = PP.period;
            do {
                if ( ++s >= PP.period )
                {
                    displayOffscreen(view, magnify);
                    if ( multi )
                        blitBuffers(fbo, multi, W, H);
                    player.saveView("movie", frm++, 2);
                    s = 0;
                }
            } while ( 0 == thread.loadNextFrame() );
        }
        
        printf("\n");
        if ( multi )
            OffScreen::releaseBuffer();
        OffScreen::releaseBuffer();
        OffScreen::closeContext();
        arg.warnings(std::cerr);
        return EXIT_SUCCESS;
    }
    
    arg.warnings(std::cerr);

    //--------- initialize Window system and create Window
#ifndef HEADLESS_PLAYER

#ifdef __APPLE__
    glutInit(&argc, argv);
#endif
    
    //register all the GLUT callback functions:
    glApp::actionFunc(processMouseClick);
    glApp::actionFunc(processMouseDrag);
    glApp::normalKeyFunc(processNormalKey);
    glApp::createWindow(displayLive);

    if ( mode == ONSCREEN_MOVIE )
    {
        PP.exit_at_eof = true;
        PP.save_images = true;
        PP.play = 1;
    }
    
    //-------- initialize graphical user interface and graphics

    try
    {
        gle::initialize();
        player.setStyle(DP.style);
        buildMenus();
        glutAttachMenu(GLUT_RIGHT_BUTTON);
        glutMenuStatusFunc(menuCallback);
        glutTimerFunc(100, timerCallback, 0);
        if ( glApp::isFullScreen() )
            glutFullScreen();
    }
    catch ( Exception & e )
    {
        print_magenta(std::cerr, "Error: "+e.brief());
        std::cerr << e.info() << '\n';
        return EXIT_FAILURE;
    }
    
    if ( player.goLive )
    {
        try
        {
            thread.period(PP.period);
            
            if ( has_frame )
                thread.extend();
            else
                thread.start();
        }
        catch( Exception & e )
        {
            print_magenta(std::cerr, "Error: "+e.brief());
            std::cerr << e.info() << '\n';
            return EXIT_FAILURE;
        }
    }

    //start the GLUT event handler:
    glutMainLoop();
#endif

    return EXIT_SUCCESS;
}
