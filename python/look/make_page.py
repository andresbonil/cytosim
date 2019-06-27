#!/usr/bin/env python
#
# make_page.py creates a HTML page with links to files in given directories
#
# Copyright F. Nedelec, 14.12.2007 -- 4.2015

"""
Synopsis:

    Generates an HTML page linking images & movies found in sub-directories

Usage:
    
    make_page.py [recursive=1] [tile=INT] [width=INT] [height=INT] [OUTPUT] DIRECTORIES

Info:
    
    If `tile` is set, a HTML table will be generated, and the value specified will be the
    number of entries in each line. The name of the ouput file ('page.html' by default), 
    can be specified on the command line (extension should be .html).
    `width` or `height` specify the size in pixels at which images will appear on the HTML.


Copyright F. Nedelec, EMBL
Created  14.12.2007
Modified 3.2010, 5.2012, 11.2012, 7.2013, 11.2013, 4.2015
"""

import sys, os, subprocess

output = 'page.html'
out    = 0
indx   = 1
imsize = ''
tile   = 0
recurs = 0


def writeHeader():
    global out;
    out.write('<html>\n')
    out.write('<head>\n')
    try:
        out.write('<title>%s</title>\n'%os.getcwd())
    except Exception:
        pass
    out.write('<meta name="ROBOTS" content="NOINDEX, NOFOLLOW">\n')
    out.write('<script type="text/JavaScript">\n')
    out.write('var zwin=0;\n')
    out.write('function mk_win()   { zwin=window.open("","see","resizable=1,width=1280,height=1024"); }\n')
    #out.write('function set(i, m)  { document.images[i].src=m; }\n')
    out.write('function zoom(m)    { mk_win(); zwin.location=m; }\n')
    #out.write('function stat(s)    { window.status=s; }\n')
    out.write('</script>\n')
    out.write('</head>\n')
    out.write('<body bgcolor="#AAAAAA">\n')
    out.write('<h2>\n')
    out.write('<a href="index.html">index</a>\n')
    out.write('</h2>\n')


def writeFooter():
    global out;
    out.write('</body>\n')
    out.write('</html>\n')


def getImageSize(file):
    """
    Call ffprobe, which is a tool that comes with ffmpeg,
    and parse output to extract size of videos or images
    """
    res = [256, 256];
    proc = subprocess.Popen(['ffprobe', '-v', 'quiet', '-show_streams', file], stdout=subprocess.PIPE)
    if not proc.wait():
        for line in proc.stdout:
            [key, equal, value] = line.partition('=')
            if key=="width":
                res[0] = int(value)
            elif key=="height":
                res[1] = int(value)
    return res


def writeImageLinks(paths):
    global out, imsize
    for path in sorted(paths):
        out.write('<a href="javascript:zoom(\'%s\');">\n' % path);
        out.write('  <img %s src="%s" alt="%s">\n' % (imsize, path, path));
        out.write('</a>\n')
    if paths:
        out.write('\n')


def writeMovieLinks(paths):
    global out
    for path in sorted(paths):
        size = getImageSize(path)
        out.write('<video controls="controls" width="%s" height="%s" loop="true" alt="%s">\n' % (size[0], size[1], path));
        out.write('  <source src="%s" type="video/mp4">\n' % path);
        out.write('  This is a HTML5 VIDEO element\n');
        out.write('</video>\n');
    if paths:
        out.write('\n')


def selectFiles(dirpath, files):
    """
    Return files to be linked in the HTML
    """
    images = []
    movies = []
    for f in files:
        path = os.path.join(dirpath, f)
        [name, ext] = os.path.splitext(f)
        if ext in ['.png', '.jpg', '.gif', '.tif', '.svg']:
            images.append(path)
        elif f.endswith('.mp4') or f.endswith('.mov'):
            movies.append(path)
    return [images, movies]


def process(dirpath, directories, files):
    """
    Write HTML code for given directory
    """
    global out, indx, tile, recurs
    dirname = dirpath.lstrip('./')
    if tile > 0:
        out.write('<td>\n')
    out.write('<h3 style="padding:3px;margin:3px"> '+dirname)
    if 'config.cym' in files:
        out.write('\n  &mdash; <a href="%s/config.cym">config</a>' % dirpath);
    if 'movie.mp4' in files:
        out.write('\n  &mdash; <a href="%s/movie.mp4">movie</a>' % dirpath);
    out.write('\n</h3>\n')
    [i, m] = selectFiles(dirpath, files)
    writeImageLinks(i)
    writeMovieLinks(m)
    if recurs:
        for d in directories:
            process_dir(os.path.join(dirpath, d))
    if tile > 0:
        out.write('</td>\n')
        if indx % tile == 0:
            out.write('</tr><tr>\n')
    indx += 1


def process_dir(dirpath):
    """call process() with appropriate arguments"""
    files = []
    directories = []
    for f in os.listdir(dirpath):
        if os.path.isdir(os.path.join(dirpath, f)):
            directories.append(f)
        else:
            files.append(f)
    process(dirpath, directories, files)

#------------------------------------------------------------------------

def main(args):
    """generates HTML page"""
    global output, out, imsize, tile, recurs
    paths = []
    
    for arg in args:
        [key, equal, value] = arg.partition('=')
        if key and equal=='=' and value:
            if key=='width' or key=='height':
                imsize = arg
            elif key=='table':
                tile = int(value)
            elif key=='tile':
                tile = int(value)
            elif key=='recursive':
                recurs = int(value)
            else:
                sys.stderr.write("ignored '%s' on command line\n" % arg)
        else:
            if os.path.isdir(arg):
                paths.append(arg)
            elif arg.endswith('.html'):
                output = arg
            else:
                sys.stderr.write("ignored '%s' on command line\n" % arg)

    if not paths:
        sys.stderr.write("You must specify a path on the command line\n")
        sys.exit()

    try:
        out = open(output, 'w')
    except Exception as e:
        sys.stderr.write("Error creating file `%s': %s\n" % (output, repr(e)));
        out = sys.stdout;
    writeHeader()
    
    if tile > 0:
        out.write('<table border="0" align="center" cellpadding="1">\n')
        out.write('<tr>\n')

    for p in paths:
        process_dir(p)
    
    if tile > 0:
        out.write('</tr>\n')
        out.write('</table>\n')
    
    writeFooter()
    if out != sys.stdout:
        out.close()
        print("generated '%s' with %i entries" % (output, indx-1))


#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv)>1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

