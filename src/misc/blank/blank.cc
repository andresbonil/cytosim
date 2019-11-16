// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/*
 This is a template for a simulation using Cytosim's libraries
*/

#include "gle.h"
#include "glut.h"
#include "glapp.h"
#include "filewrapper.h"
#include "blank_param.h"
#include "glossary.h"
#include "random.h"

#include "thing.h"

// parameter list:
BlankParam PAM;

// random number generator

//------------------------------------------------------------------------------
#include <pthread.h>

pthread_t       slave;
pthread_mutex_t mutex;
pthread_cond_t  condition;

// simulation time:
real sTime = 0;

// objects:
const unsigned int maxThings = 1024;
Thing things[maxThings];


//------------------------------------------------------------------------------
void simReset()
{
    RNG.seed();

    if ( PAM.max > maxThings ) 
        PAM.max = maxThings;
    
    sTime = 0;
    for(unsigned i = 0; i < PAM.max; ++i)
        things[i].reset();
}

void simStep()
{
    sTime += PAM.time_step;
    for(unsigned i = 0; i < PAM.max; ++i)
        things[i].step();
    glApp::postRedisplay();
}


void * simulateForever(void *)
{
    simReset();
    while ( 1 )
    {
        for (unsigned r = 0; r < PAM.repeat; ++r )
            simStep();
        pthread_cond_wait(&condition, &mutex);
    }
}


void simulateToFile(FILE * out)
{
    for(unsigned i = 0; i < PAM.max; ++i)
        things[i].write(out, sTime);

    for(unsigned f = 1; f <= PAM.repeat; ++f )
    {
        while ( sTime < f * PAM.time )
            simStep();
        for (unsigned i = 0; i < PAM.max; ++i )
            things[i].write(out, sTime);
        fprintf(out, "\n");
    }
}

//------------------------------------------------------------------------------
void processNormalKey(unsigned char c, int x, int y)
{
    switch (c)
    {
        case 'p':
        {
            const unsigned int minDelay = 2;
            if ( PAM.delay < 2*minDelay )
            {
                PAM.delay = minDelay;
                PAM.repeat *= 2;
                if ( PAM.repeat > 1024 )
                    PAM.repeat = 1024;
            }
            else {
                PAM.delay /= 2;
            }
            glApp::flashText("Delay %i ms, %i steps/display", PAM.delay, PAM.repeat);
        } break;
        case 'o':
            if ( PAM.repeat > 1 ) {
                PAM.repeat /= 2;
            }
            else {
                PAM.repeat = 1;
                PAM.delay = ( PAM.delay < 1024 ) ? 2*PAM.delay : 2048;
            }
            glApp::flashText("Delay %i ms, %i steps/display", PAM.delay, PAM.repeat);
            break;
        case 's':
            PAM.repeat = !PAM.repeat;
            break;
        case 'Z':
            simReset();
            glApp::flashText("Reset simulation");
            break;
        case 'z':
            glApp::resetView();
            glApp::flashText("Reset view");
            break;
        case 'r':
            for(unsigned i = 0; i < PAM.max; ++i)
                things[i].write(stdout, sTime);
            break;
        case 'R':
            PAM.write(std::cout);
            break;
        default:
            glApp::processNormalKey(c,x,y);
            return;
    }

    glutPostRedisplay();
}

/// mouse callback for shift-click, with unprojected mouse position
void processMouseClick(int, int, const Vector3& a, int)
{
    
    glApp::flashText("click %.1f %.1f %.1f", a.XX, a.YY, a.ZZ);
    glApp::postRedisplay();
}

/// mouse callback for shift-drag, with unprojected mouse positions
void processMouseDrag(int, int, Vector3& a, const Vector3& b, int)
{
    
    glApp::flashText("move %.1f %.1f %.1f", b.XX, b.YY, b.ZZ);
    glApp::postRedisplay();
}


///timer callback
void timerFunction(int value)
{
    glutTimerFunc(PAM.delay, timerFunction, 1);
    pthread_cond_signal(&condition);
}

//------------------------------------------------------------------------------
void display(View& view, int)
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    
    view.setLabel(std::to_string(sTime));

    glPointSize(6);
    glBegin(GL_POINTS);
    pthread_mutex_lock(&mutex);
    for(unsigned i = 0; i < PAM.max; ++i)
        things[i].draw();
    pthread_mutex_unlock(&mutex);
    glEnd();
}

//------------------------------------------------------------------------------
void help()
{
    printf("Blank 1.0 by Francois Nedelec, Copyright EMBL 2007\n"
           "www.cytosim.org\n\n"
           "Command-line options:\n"
           "parameters  display list of parameters\n"
           "help        display this help\n"
           "keys        display list of keyboard controls\n"
           "P=###       set value of parameter `P'\n");
}

void help_keys(std::ostream& os = std::cout)
{
    os << "Blanksim 0.1,  Francois J. Nedelec 2007-2011\n\n";
    glApp::help(os);
    os << "\n Keyboard Controls:\n";
    os << " p o        Increase/decrease display speed\n";
    os << " s          Stop/start simulation\n";
    os << " f          toggle full-screen display\n";
    os << " x          Display XYZ axis\n";
    os << " Z z        Reset simulation, reset view\n";
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    //parse the config file
    Glossary arg;
    if ( arg.read_strings(argc-1, argv+1) )
        return EXIT_FAILURE;
    PAM.read(arg);
    
    if ( arg.use_key("help") )
    {
        help();
        help_keys();
        return EXIT_SUCCESS;
    }
    
    if ( arg.use_key("parameters") )
    {
        PAM.write(std::cout);
        return EXIT_SUCCESS;
    }
    
    glutInit(&argc, argv);
    glApp::setDimensionality(2);
    glApp::setScale(10);
    glApp::attachMenu(GLUT_RIGHT_BUTTON);
    glApp::normalKeyFunc(processNormalKey);
    glApp::actionFunc(processMouseClick);
    glApp::actionFunc(processMouseDrag);

    if ( PAM.config.length() )
        arg.read_file(PAM.config);
    
    PAM.read(arg);
    glApp::currentView().read(arg);
    arg.warnings(std::clog);
    
    glApp::createWindow(display);
    gle::initialize();
    simReset();
    
    pthread_mutex_init(&mutex, nullptr);
    pthread_cond_init(&condition, nullptr);
    pthread_create(&slave, nullptr, &simulateForever, nullptr);

    timerFunction(0);
    glutMainLoop();
    return EXIT_SUCCESS;
}
