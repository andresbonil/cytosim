// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SAVEIMAGE_H
#define SAVEIMAGE_H

#include <cstdio>


#ifndef GL_VERSION_1_1

typedef unsigned char GLubyte;

#endif


/// Can save pixel array to files in PNG or PPM format
/**
 - PPM files do not require any library, and thus writing is supported on any platform.
   They can be read by various software, in particular the 'pbmplus' toolkit
   However, they are big and usually not supported by most viewers.
   Known viewers on Mac OS X: ImageJ, GraphicConverter, ToyViewer.
 - PNG is a modern all-purpose file format that supports RGBA format, and is thus
   very well suited to export OpenGL scenes.
 .
 
 PNG support requires 'libpng', which must be installed separately.
 */
namespace SaveImage
{
    /// destination of error messages (set to zero to suppress output)
    static FILE * err = stderr;
    
    /// open a file for binary write
    FILE * openFile(const char * name);

    /// Netpbm pixel image format, 3 one-byte componentss per pixels (R, G, B)
    int saveColorPPM(FILE*, const GLubyte[], int width, int height);
    
    /// write PNG image
    int savePNG(FILE*, const GLubyte[], int bit_depth, int nb_colors, int width, int height);

    /// Portable Network Graphic format, 4 one-byte components per pixels (R, G, B, A)
    int saveAlphaPNG(FILE*, const GLubyte[], int width, int height);

    /// Portable Network Graphic format, 3 one-byte components per pixels (R, G, B)
    int saveColorPNG(FILE*, const GLubyte[], int width, int height);

    /// gray-level PNG format with one 2-byte component per pixels
    int saveGrayPNG(FILE*, const GLubyte[], int width, int height);
    
    /// save pixels[] and return error-code
    int savePixels(FILE*, const char format[], const GLubyte[], int width, int height);
    
    /// save pixels[] and return error-code
    int savePixels(FILE*, const char format[], const GLubyte[], int width, int height, int downsample);

    /// save pixels[] and return error-code
    int savePixels(const char* name, const char format[], const GLubyte[], int width, int height, int downsample);
    
    /// downsample `src` to set `dst`
    void downsampleRGB(GLubyte* dst, GLubyte const* src, unsigned W, unsigned H, unsigned bin);

    //-------------------- use the functions below: ---------------------
    
    /// true if 'format' is the 3-letter file-entension of a supported image format
    /**
     'ppm' is always supported, and 'png' may be supported depending on compilation
     The 3-letter format string can be lowercase or uppercase.
     */
    bool supported(const char format[]);

    /// save entire viewport in a new file called 'name'. Returns error-code
    int readPixels(int x, int y, int width, int height, void * pixels);

    /// save entire viewport in a new file called 'name'. Returns error-code
    int saveEntireImage(const char* name, const char format[], int downsample=1);

    /// save a region of the current buffer in a new file called 'name'. Returns error-code
    int saveImage(const char* name, const char format[], const int viewport[], int downsample=1);

     /// save an image with higher resolution (better version)
    int saveMagnifiedImage(int mag, const char* name, const char format[], int width, int height, void (*display)(int, void *), void * arg, int downsample=1);

    /// save an image with higher resolution
    int saveCompositeImage(int mag, const char* name, const char format[], int width, int height, double pixel_size, void (*display)(int, void *), void * arg, int downsample=1);
}


#endif
