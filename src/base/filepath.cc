// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "filepath.h"
#include <sys/param.h>
#include <sys/stat.h>
#include <cstdlib>
#include <libgen.h>
#include <unistd.h>
#include <dirent.h>
#include "exceptions.h"


std::string FilePath::get_cwd()
{
    char cwd[1024] = { 0 };
    if ( getcwd(cwd, sizeof(cwd)) )
        return cwd;
    else
        return "";
}


bool FilePath::is_file(const char path[])
{
    struct stat s;
    if ( 0 == stat(path, &s) )
        return S_ISREG(s.st_mode);
    return false;
}


bool FilePath::is_dir(const char path[])
{
    struct stat s;
    if ( 0 == stat(path, &s) )
        return S_ISDIR(s.st_mode);
    return false;
}


int FilePath::make_dir(const char name[])
{
    return mkdir(name, S_IRWXU|S_IRWXG|S_IRWXG|S_IXOTH);
}


int FilePath::change_dir(const char path[], bool make)
{
    if ( !path || *path == 0 )
        return 0;
    if ( make )
        mkdir(path, S_IRWXU|S_IRWXG|S_IRWXG|S_IXOTH);
    return chdir(path);
}


std::vector<std::string> FilePath::list_dir(const char path[])
{
    std::vector<std::string> res;
    DIR * dp = opendir(path);
    if ( dp )
    {
        struct dirent * ep = readdir(dp);
        while ( ep )
        {
            res.push_back(ep->d_name);
            ep = readdir(dp);
        }
        
        closedir(dp);
    }
    return res;
}


std::vector<std::string> FilePath::list_dir(const char path[], std::string const& ext)
{
    std::vector<std::string> res;
    DIR * dp = opendir(path);
    if ( dp )
    {
        struct dirent * ep = readdir(dp);
        while ( ep )
        {
            std::string name(ep->d_name);
            if ( name.size() > ext.size()
                &&  0==name.compare(name.size()-ext.size(), ext.size(), ext) )
                res.push_back(ep->d_name);
            ep = readdir(dp);
        }
        
        closedir(dp);
    }
    
#if ( 0 )
    for ( size_t i = 0; i < res.size(); ++i )
        std::clog << "   " << res[i] << std::endl;
#endif
    return res;
}


std::string FilePath::dir_part(const char path[])
{
    char* res, tmp[PATH_MAX] = { 0 };
    
    if ( !realpath(path, tmp) )
        return ".";
    
    res = dirname(tmp);
    
    if ( !res )
        throw InvalidIO("FilePath: stdlib::dirname() failed");
    
    return res;
}


std::string FilePath::file_part(const char path[])
{
    char str[MAXPATHLEN] = { 0 };
    strncpy(str, path, MAXPATHLEN);
    char * res = basename(str);
    
    if ( !res )
        throw InvalidIO("FilePath: stdlib::basename() failed");
    
    return std::string(res);
}


std::string FilePath::full_name(std::string const& dir, std::string const& file)
{
    //if a full path is already specified, we do nothing
    if ( dir.size() > 0  &&  file.size() > 0  &&  file[0] != '/' )
    {
        std::string res = dir + file;
        
        //remove trailling '/' if present
        size_t s = res.size();
        if ( 1 <= s  &&  res[s-1] == '/' )
            res[s-1] = '\0';
        
        return res;
    }
    return file;
}

