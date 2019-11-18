// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <cctype>
#include "point_disp.h"
#include "offscreen.h"
#include "glossary.h"
#include "opengl.h"
#include "glut.h"
#include "gle.h"


/// if this is defined, the pixelmap are stored in graphical memory
#define POINTDISP_USES_PIXEL_BUFFERS 0


void PointDisp::clearPixelmaps()
{
#if POINTDISP_USES_PIXELMAPS
    nPix = 0;
    for ( int i = 0; i < 3; ++i )
    {
        bmp[i] = nullptr;
        pbo[i] = 0;
    }
#endif
}


PointDisp::PointDisp(const std::string& k, const std::string& n)
: Property(n)
{
    mKind = k;
    clearPixelmaps();
    clear();
}


PointDisp::PointDisp(PointDisp const& o) : Property(o)
{
    mKind        = o.mKind;
    visible      = o.visible;
    color        = o.color;
    color2       = o.color2;
    coloring     = o.coloring;
    size         = o.size;
    width        = o.width;
    shape        = o.shape;
    style        = o.style;
    symbol       = o.symbol;
    symbol_color = o.symbol_color;

    clearPixelmaps();
}


PointDisp& PointDisp::operator = (PointDisp const& o)
{
    mKind        = o.mKind;
    visible      = o.visible;
    color        = o.color;
    color2       = o.color2;
    coloring     = o.coloring;
    size         = o.size;
    width        = o.width;
    shape        = o.shape;
    style        = o.style;
    symbol       = o.symbol;
    symbol_color = o.symbol_color;
    
    clearPixelmaps();
    
    return *this;
}


PointDisp::~PointDisp()
{
#if POINTDISP_USES_PIXELMAPS
    releasePixelmap();
    
#if POINTDISP_USES_PIXEL_BUFFERS
    if ( pbo[0] > 0 )
    {
        glDeleteBuffers(3, pbo);
        pbo[0] = 0;
        pbo[1] = 0;
        pbo[2] = 0;
    }
#endif
#endif
}


void PointDisp::clear()
{
    visible      = 1;
    color        = 0x888888FF;
    color2       = 0x777777FF;
    coloring     = 0;
    size         = 5;
    width        = 2;
    shape        = 'o';
    style        = 7;
    symbol       = 0;
    symbol_color = 0xFFFFFFFF;
}


void PointDisp::strokeShape() const
{
    switch ( tolower(shape) )
    {
        case 'v': gle::gleNablaL();     break;
        case 't': gle::gleTriangleL();  break;
        case 'q': gle::gleSquareL();    break;
        case 'r': gle::gleRectangleL(); break;
        case 'p': gle::glePentagonL();  break;
        case 'h': gle::gleHexagonL();   break;
        case 's': gle::gleStarL();      break;
        case '+': gle::glePlusL();      break;
        case 'c': gle::gleCircle();    break;
        default: break;
    }
}


void PointDisp::paintShape() const
{
    switch ( tolower(shape) )
    {
        case 'v': gle::gleNablaS();     break;
        case 't': gle::gleTriangleS();  break;
        case 'q': gle::gleSquareS();    break;
        case 'r': gle::gleRectangleS(); break;
        case 'p': gle::glePentagonS();  break;
        case 'h': gle::gleHexagonS();   break;
        case 's': gle::gleStarS();      break;
        case '+': gle::glePlusS();      break;
        case 'c': break;
        default:  gle::gleDisc();    break;
    }
}

void PointDisp::strokeA() const
{
    paintShape();
    if ( width > 0.5 )
    {
        //draw a bright rim
        color.darken(2.0).load();
        glLineWidth(3.0);
        strokeShape();
    }
    
    if ( symbol )
    {
        glScalef(1.0f/80, 1.0f/80, 1);
        /*  glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, C) stokes
         character C of width ~104.76 units, and ~150 unit high max
         The translation brings it near the center. */
        if ( islower(symbol) )
            glTranslatef(-52.35f, -35, 0);
        else
            glTranslatef(-52.35f, -50, 0);
        symbol_color.load();
        glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, symbol);
    }
}

void PointDisp::strokeI() const
{
    paintShape();
    
    // radius of the spot that indicate an inactive Hand
    const GLfloat DOT_SIZE = 0.55f;

    // draw a transparent hole in the center:
    glPushAttrib(GL_COLOR_BUFFER_BIT|GL_ENABLE_BIT);
    glPushMatrix();
    glScalef(DOT_SIZE, DOT_SIZE, DOT_SIZE);
    glDisable(GL_BLEND);
    glDisable(GL_ALPHA_TEST);
    glColor4f(0, 0, 0, 0);
    gle::gleDiscB(); //strokeShape();
    glPopMatrix();
    glPopAttrib();
}


#pragma mark - Bitmaps

#if POINTDISP_USES_PIXELMAPS


/**
 Allocate memory for pixelmaps of size `pixSize x pixSize`
 3 bitmaps x 4 colors x pixSize x pixSize pixels
 \todo use Graphic card memory for the bitmaps
 */
void PointDisp::allocatePixelmap()
{
    mOffs = -0.5f * pixSize;
    
    // allocate only if needed
    if ( !bmp[0] || pixSize != nPix )
    {
        if ( bmp[0] )
            delete(bmp[0]);
    
        unsigned dd = pixSize * pixSize;
        GLubyte * mem = new GLubyte[12*dd];
        bmp[0] = mem;
        bmp[1] = mem + 4*dd;
        bmp[2] = mem + 8*dd;
        
        for ( unsigned y = 0; y < 12*dd; ++y )
            mem[y] = 0;
        
        nPix = pixSize;
    }
}


void PointDisp::releasePixelmap()
{
    if ( bmp[0] )
    {
        delete(bmp[0]);
        bmp[0] = nullptr;
    }
    nPix = 0;
}


/// print pixel map in ASCII
void printPixels(GLubyte const* pix, unsigned sx, unsigned sy)
{
    for ( unsigned y = 0; y < sy; ++y )
    {
        for ( unsigned x = 0; x < sx; ++x )
            printf("%02X", pix[4*(x+sx*y)+3]);
        printf("\n");
    }
}

/**
 This will downsample pixelmap `src` and set destination `dst`. The pixel
 array `dst` should be of size `4*sx*sy` with 4 bytes per pixels: R, G, B and A,
 while `src` should be `bin*bin` times larger. Pixels are stored in row order
 from the lowest to the highest row, left to right in each row (as in OpenGL).
 The pixels components of `src` are averaged to produce `dst`.
 */
void PointDisp::downsampleRGBA(GLubyte dst[], unsigned sx, unsigned sy,
                               GLubyte const src[], unsigned bin)
{
    const unsigned s = bin * bin;

#if ( 0 )
    //reset destination:
    for ( unsigned u = 0; u < sx*sy; ++u )
    {
        dst[4*u  ] = 0xFF;
        dst[4*u+1] = 0xFF;
        dst[4*u+2] = 0xFF;
        dst[4*u+3] = 0xFF;
    }
#endif
    
    for ( unsigned x = 0; x < sx; ++x )
    for ( unsigned y = 0; y < sy; ++y )
    {
        unsigned r = 0, g = 0, b = 0, a = 0;
        for ( unsigned dx = 0; dx < bin; ++dx )
        for ( unsigned dy = 0; dy < bin; ++dy )
        {
            GLubyte const* p = src + 4 * ( dx+bin*(x+sx*(dy+bin*y)) );
            r += p[0];
            g += p[1];
            b += p[2];
            a += p[3];
        }
            
        dst[4*(x+sx*y)  ] = (GLubyte)( r / s );
        dst[4*(x+sx*y)+1] = (GLubyte)( g / s );
        dst[4*(x+sx*y)+2] = (GLubyte)( b / s );
        dst[4*(x+sx*y)+3] = (GLubyte)( a / s );
    }
}


void PointDisp::storePixelmap(GLubyte* bitmap, unsigned dim, GLuint pbi) const
{
#if POINTDISP_USES_PIXEL_BUFFERS
    assert_true(pbi);
    glBindBuffer(GL_PIXEL_PACK_BUFFER_ARB, pbi);
    glBufferData(GL_PIXEL_PACK_BUFFER_ARB, 4*dim*dim, bitmap, GL_STATIC_DRAW);
    glBindBuffer(GL_PIXEL_PACK_BUFFER_ARB, 0);
    assert_true(glIsBuffer(pbi));
    //gle::gleReportErrors(stderr, "PointDisp::storePixelmap");
#endif
}


#if ( 0 )

#include "saveimage.h"

// Export bitmap to file in PNG format
void PointDisp::savePixelmap(GLubyte* bitmap, unsigned dim, unsigned id) const
{
    if ( SaveImage::supported("png") )return given name of property
    {
       char str[32];
        snprintf(str, sizeof(str), "bitmap_%s_%02u.png", name_str(), id);
        FILE * file = fopen(str, "w");
        if ( file && !ferror(file ) )
        {
            SaveImage::saveAlphaPNG(file, bitmap, dim, dim);
            fclose(file);
            std::clog << "PointDisp saved " << str << '\n';
        }
    }
}
#endif


void PointDisp::drawPixelmap(Vector const& pos, unsigned ii) const
{
    gle::gleRasterPos(pos);
    //translate to center the bitmap:
    glBitmap(0,0,0,0,mOffs,mOffs,nullptr);
#if POINTDISP_USES_PIXEL_BUFFERS
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, pbo[ii]);
    glDrawPixels(nPix, nPix, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
#else
    glDrawPixels(nPix, nPix, GL_RGBA, GL_UNSIGNED_BYTE, bmp[ii]);
#endif
    gle::gleReportErrors(stderr, "PointDisp::drawPixelmap");
}


/**
 `sampling` defines the level of oversampling used to improve the quality of bitmaps
 */
void PointDisp::makePixelmaps(GLfloat uFactor, unsigned sampling)
{
    assert_true(pixSize==nPix);
    //gle::gleReportErrors(stderr, "1 PointDisp::makePixelmaps");
    
    glPushAttrib(GL_PIXEL_MODE_BIT|GL_VIEWPORT_BIT|GL_ENABLE_BIT|GL_COLOR_BUFFER_BIT);

    unsigned dim = sampling * nPix;
#ifdef __APPLE__
    GLuint offscreen = OffScreen::createBuffer(dim, dim, 0);
#else
    GLuint offscreen = 0;
#endif
    GLfloat s = 0.5f * sampling * size * uFactor;
    GLfloat t = 0.5f * dim;
    GLfloat w = width * sampling * uFactor;
    GLint vp[4];

    glDisable(GL_LIGHTING);
    glDisable(GL_BLEND);
    glDisable(GL_ALPHA_TEST);
    
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    
    //match projection to viewport:
    glGetIntegerv(GL_VIEWPORT, vp);
    glOrtho(0, vp[2], 0, vp[3], 0, 1);
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    glTranslatef(t, t, 0);
    glScalef(s,s,s);
    if ( w > 0.0 ) glLineWidth(w);
    // we use a transparent background, because points will overlap
    glClearColor(0,0,0,0);
    
    for ( int i = 0; i < 3; ++i )
    {
        glClear(GL_COLOR_BUFFER_BIT);
        switch ( i )
        {
            case 0:
                color2.load();
                strokeI();
                break;
            case 1:
                color2.load();
                strokeA();
                break;
            case 2:
                color.load();
                strokeA();
                break;
        }
        if ( sampling > 1 )
        {
            GLubyte * tmp = new GLubyte[4*dim*dim];
            glReadPixels(0, 0, dim, dim, GL_RGBA, GL_UNSIGNED_BYTE, tmp);
            downsampleRGBA(bmp[i], nPix, nPix, tmp, sampling);
            delete[] tmp;
#if ( 0 )
            //savePixelmap(tmp, dim, i);
            std::clog << name() << i << "\n";
            printPixels(bmp[i], nPix, nPix);
#endif
        }
        else
        {
            glReadPixels(0, 0, nPix, nPix, GL_RGBA, GL_UNSIGNED_BYTE, bmp[i]);
            //savePixelmap(bmp[i], nPix, i+10);
        }
        //gle::gleReportErrors(stderr, "5 PointDisp::makePixelmaps");
        storePixelmap(bmp[i], nPix, pbo[i]);
    }
 
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    
    // release offscreen if used
    if ( offscreen )
        OffScreen::releaseBuffer();

    glPopAttrib();
}

#else

/// draw inactive state
void PointDisp::drawI(Vector const& pos) const
{
    if ( perceptible )
    {
        glPushMatrix();
        gle::gleTranslate(pos);
        gle::gleScale(realSize);
        color2.load();
        gle::gleDisc();
        glPopMatrix();
    }
}

/// draw active state, unattached
void PointDisp::drawF(Vector const& pos) const
{
    if ( perceptible )
    {
        glPushMatrix();
        gle::gleTranslate(pos);
        gle::gleScale(realSize);
        color2.load();
        strokeA();
        glPopMatrix();
    }
}

/// draw active state, attached
void PointDisp::drawA(Vector const& pos) const
{
    if ( perceptible )
    {
        glPushMatrix();
        gle::gleTranslate(pos);
        gle::gleScale(realSize);
        color.load();
        strokeA();
        glPopMatrix();
    }
}

#endif


void PointDisp::prepare(GLfloat uf, GLfloat sf, bool make_maps)
{
    realSize    = size * sf;
    unsigned sz = (unsigned)ceil(uf*(size+width));
    // use a multiple of 4 pixels:
    pixSize     = ( sz + 4 ) & ~3;
    perceptible = visible && ( uf*(size+width) > 0.25 );
    
    if ( make_maps )
    {
#if POINTDISP_USES_PIXELMAPS
#if POINTDISP_USES_PIXEL_BUFFERS
        if ( pbo[0] == 0 )
            glGenBuffers(3, pbo);
#endif
        
        if ( pixSize != nPix )
        {
            //gle::gleReportErrors(stderr, "1 PointDisp::prepare");
            allocatePixelmap();
            //fprintf(stderr, " new %i bitmap for %s\n", pixSize, name_str());
            //gle::gleReportErrors(stderr, "2 PointDisp::prepare");
        }
        
        makePixelmaps(uf, 3);
#endif
    }
}


#pragma mark - I/O


void PointDisp::read(Glossary& glos)
{
    glos.set(visible,      "visible");
    
    // set 'color2' as a darker tone of 'color':
    if ( glos.set(color,   "color") )
        color2 = color.alpha(0.5);
    glos.set(color2,       "color", 1) || glos.set(color2, "back_color");
    glos.set(coloring,     "coloring");
    
    // if 'size' is specified, width is set accordingly:
    if ( glos.set(size,    "size") )
        width = 2 * size / 3;

    // alternative syntax:
    glos.set(size,         "point_size");
#ifdef BACKWARD_COMPATIBILITY
    glos.set(size,         "points");
    glos.set(shape,        "points", 1);
#endif

    glos.set(width,        "width");
    glos.set(style,        "style");
    glos.set(shape,        "shape");
    glos.set(symbol,       "symbol");
    glos.set(symbol_color, "symbol", 1);
    
    if ( ! isprint(symbol) )
        symbol = 0;
    
#if POINTDISP_USES_PIXELMAPS
    releasePixelmap();
#endif
}


void PointDisp::write_values(std::ostream& os) const
{
    write_value(os, "visible",     visible);
    if ( color2 != color.alpha(0.5) )
        write_value(os, "color",   color, color2);
    else
        write_value(os, "color",   color);
    write_value(os, "coloring",    coloring);
    write_value(os, "size",        size);
    write_value(os, "width",       width);
    write_value(os, "shape",       shape);
    write_value(os, "style",       style);
    write_value(os, "symbol",      symbol, symbol_color);
}

