// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/**
 'Frametool' is a simple utility that can read and extract frames in
 Cytosim's trajectory files, usually called "objects.cmo".
 
 It only uses the START and END tags of frames, and does not care
 about the organization of the data contained between these tags.
 
 'Frametool' can extract frames from the file, which is useful
 for example to reduce the size of 'objects.cmo' by dropping some frames.

 You can reduce the file size by half by dropping every odd frame:
 > frametool objects.cmo 0:2: > o.cmo
 > mv o.cmo objects.cmo
 
 Another tool 'sieve' can be used to read/write object-files,
 allowing finer manipulation of the simulation frames.
*/

#include <errno.h>
#include <cstdio>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>


enum { COUNT, COPY, LAST, SIZE, EPID, SPLIT };
enum { UNKNOWN, FRAME_START, FRAME_SECTION, FRAME_END };

FILE * output = stdout;
const size_t buf_size = 64;
char buf[buf_size];

unsigned long frame_pid = 0;


FILE * openfile(char name[], char const* mode)
{
    FILE * file = fopen(name, mode);
    
    if ( file==nullptr )
    {
        fprintf(stderr, "Could not open file `%s'\n", name);
        return nullptr;
    }
    
    if ( ferror(file) )
    {
        fclose(file);
        fprintf(stderr, "Error opening file `%s'\n", name);
        return nullptr;
    }
    
    return file;
}


/**
 read a line, and returns a code indicating if this is the start
 or the end of a cytosim frame
 */
int whatline(FILE* in, FILE* out)
{
    char *const end = buf + buf_size - 1;
    char * ptr = buf;

    int c = 0;
    do {
        c = getc_unlocked(in);
        
        if ( c == EOF )
            return EOF;

        if ( ptr < end )
            *ptr++ = (char)c;
        
        if ( out )
            putc_unlocked(c, out);
        
    } while ( c != '\n' );
    
    // fill-in with zeros:
    if ( ptr < end )
        *ptr = 0;
    else
        *end = 0;
    
    if ( *buf == '#' )
    {
        if ( 0 == strncmp(buf, "#frm ", 5) )     return FRAME_START;
        if ( 0 == strncmp(buf, "#frame ", 7) )   return FRAME_START;
        if ( 0 == strncmp(buf, "#Cytosim ", 9) )
        {
            frame_pid = strtoul(buf+10, nullptr, 10);
            return FRAME_START;
        }
        if ( 0 == strncmp(buf, "#end ", 5) )     return FRAME_END;
        if ( 0 == strncmp(buf, " #end ", 6) )    return FRAME_END;
        if ( 0 == strncmp(buf, "#section ", 9) ) return FRAME_SECTION;
    }
    return UNKNOWN;
}


//=============================================================================

void error(const char* message)
{
    fprintf(stderr, "ERROR: %s\n", message);
    exit(EXIT_FAILURE);
}


/// Slice represents a regular subset of indices
class Slice
{
    unsigned s; ///< start
    unsigned i; ///< increment
    unsigned e; ///< end
    
public:
    
    Slice()
    {
        s =  0;
        i =  1;
        e = ~0U;
    }
    
    Slice(const char arg[])
    {
        s =  0;
        i =  1;
        e = ~0U;

        int c = 0;
        c = sscanf(arg, "%u:%u:%u", &s, &i, &e);
        //fprintf(stderr, "%s:%i\n", arg, c);
        
        if ( arg[strlen(arg)-1] == ':' )
        {
            if ( c == 3 )
                error("unexpected third ':'");
        }
        else
        {
            if ( c == 1 )
                e = s;
            if ( c == 2 )
            { e = i; i = 1; }
        }
        //fprintf(stderr, "slice %u:%u:%u\n", s, p, e);
    }
    
    bool match(unsigned n)
    {
        if ( n < s )
            return false;
        if ( e < n )
            return false;
        return 0 == ( n - s ) % i;
    }
    
    unsigned last()
    {
        return e;
    }
};

//=============================================================================

void countFrame(const char str[], FILE* in)
{
    int  frm = 0;
    int code = 0;
    do {
        code = whatline(in, nullptr);
        if ( code == FRAME_END )
            ++frm;
    } while ( code != EOF );
    
    printf("%40s: %i frames\n", str, frm);
}


void sizeFrame(FILE* in)
{
    int  code = 0, frm = -1, cnt = 0, oldcnt = 0;

    while ( code != EOF )
    {
        ++cnt;
        code = whatline(in, nullptr);
        
        if ( code == FRAME_END )
        {
            printf("%lu  frame %5i: %7i lines (%+i)\n", frame_pid, frm, cnt, cnt-oldcnt);
            oldcnt = cnt;
        }
        
        if ( code == FRAME_START )
        {
            ++frm;
            cnt = 0;
        }
    }
}


void extract(FILE* in, FILE* out, Slice sli)
{
    unsigned frm = 0;
    int  code = 0;
    FILE * file = sli.match(0) ? out : nullptr;

    while ( code != EOF )
    {
        code = whatline(in, file);
        
        if ( code == FRAME_START )
        {
            if ( frm > sli.last() )
                return;
            file = sli.match(frm) ? out : nullptr;
        }

        if ( code == FRAME_END )
        {
            if ( ++frm > sli.last() )
                return;
            file = sli.match(frm) ? out : nullptr;
        }
    }
}


void extract_pid(FILE* in, unsigned long pid)
{
    int code = 0;
    FILE * out = nullptr;

    while ( code != EOF )
    {
        code = whatline(in, out);
        
        if ( code == FRAME_START && pid == frame_pid )
            out = stdout;
        else
            out = nullptr;
    }
}


void extractLast(FILE* in)
{
    fpos_t pos, start;
    fgetpos(in, &start);

    int code = 0;
    while ( code != EOF )
    {
        code = whatline(in, nullptr);
        if ( code == FRAME_END )
        {
            start = pos;
            fgetpos(in, &pos);
        }
    }
    
    clearerr(in);
    fsetpos(in, &start);
    
    int c = 0;
    while ( 1 )
    {
        c = getc_unlocked(in);
        if ( c == EOF )
            break;
        putchar(c);
    }
    
    putchar('\n');
}


void split(FILE * in)
{
    int frm = 0;
    int  code = 0;
    char name[128] = { 0 };
    snprintf(name, sizeof(name), "objects%04i.cmo", frm);
    FILE * out = fopen(name, "w");
   
    while ( code != EOF )
    {
        code = whatline(in, out);
        
        if ( code == FRAME_END )
        {
            if ( out )
            {
                funlockfile(out);
                fclose(out);
            }
            ++frm;
            snprintf(name, sizeof(name), "objects%04i.cmo", frm);
            out = openfile(name, "w");
            if ( out )
                flockfile(out);
        }
    }
}

//=============================================================================

void help()
{
    printf("Synopsis:\n");
    printf("    `frametool` can list the frames present in a trajectory file,\n");
    printf("     and extract selected ones\n");
    printf("Syntax:\n");
    printf("    frametool FILENAME \n");
    printf("    frametool FILENAME size\n");
    printf("    frametool FILENAME pid=PID\n");
    printf("    frametool FILENAME INDICES\n");
    printf("    frametool FILENAME split\n");
    printf(" where INDICES specifies an integer or a range of integers as:\n");
    printf("        INDEX\n");
    printf("        START:END\n");
    printf("        START:\n");
    printf("        START:INCREMENT:END\n");
    printf("        START:INCREMENT:\n");
    printf("        last\n");
    printf(" The option 'split' will create one file for each frame in the input\n");
    printf("Examples:\n");
    printf("    frametool objects.cmo 0:2:\n");
    printf("    frametool objects.cmo 0:10\n");
    printf("    frametool objects.cmo last\n");
    printf("    frametool objects.cmo split\n");
}


bool is_file(const char path[])
{
    struct stat s;
    if ( 0 == stat(path, &s) )
        return S_ISREG(s.st_mode);
    return false;
}

bool is_dir(const char path[])
{
    struct stat s;
    if ( 0 == stat(path, &s) )
        return S_ISDIR(s.st_mode);
    return false;
}

int main(int argc, char* argv[])
{
    int mode = COUNT;
    char cmd[256] = "";
    char filename[256] = "objects.cmo";
    unsigned long pid = 0;

    for ( int i = 1; i < argc ; ++i )
    {
        if ( 0 == strncmp(argv[i], "help", 4) )
        {
            help();
            return EXIT_SUCCESS;
        }
        if ( is_file(argv[i]) )
            strncpy(filename, argv[i], sizeof(filename));
        else if ( is_dir(argv[i]) )
            snprintf(filename, sizeof(filename), "%s/objects.cmo", argv[i]);
        else
        {
            strncpy(cmd, argv[i], sizeof(cmd));
            
            if ( isdigit(*cmd) )
                mode = COPY;
            else if ( 0 == strncmp(cmd, "last", 4) )
                mode = LAST;
            else if ( 0 == strncmp(cmd, "split", 5) )
                mode = SPLIT;
            else if ( 0 == strncmp(cmd, "size", 4) )
                mode = SIZE;
            else if ( 0 == strncmp(cmd, "count", 5) )
                mode = COUNT;
            else if ( 0 == strncmp(cmd, "pid=", 4) )
            {
                mode = EPID;
                errno = 0;
                pid = strtoul(cmd+4, nullptr, 10);
                if ( errno )
                {
                    fprintf(stderr, "unexpected syntax");
                    return EXIT_FAILURE;
                }
            }
            else
            {
                fprintf(stderr, "unexpected command (for help, invoke `frametool help`)\n");
                return EXIT_FAILURE;
            }
        }
    }
    
    //----------------------------------------------
    
    FILE * file = openfile(filename, "r");
    
    if ( !file )
        return EXIT_FAILURE;
    
    flockfile(file);
    switch(mode)
    {
        case COUNT:
            countFrame(filename,file);
            break;

        case SIZE:
            sizeFrame(file);
            break;

        case COPY:
            extract(file, stdout, Slice(cmd));
            break;
            
        case LAST:
            extractLast(file);
            break;
            
        case EPID:
            extract_pid(file, pid);
            break;
            
        case SPLIT:
            split(file);
            break;
    }
    funlockfile(file);
    fclose(file);
    
    return EXIT_SUCCESS;
}