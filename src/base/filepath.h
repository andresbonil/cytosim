// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include <vector>
#include <string>

/// Functions to handle UNIX-like file/directory paths
namespace FilePath
{
    /// return path to current working directory
    std::string get_cwd();

    /// true if 'path' is an existing readible file
    bool is_file(const char path[]);
    
    /// true if 'path' is an existing directory
    bool is_dir(const char path[]);
    
    /// change current working directory
    int  change_dir(const char path[]);

    /// true if 'path' is an existing directory
    inline bool is_dir(std::string const& str) { return is_file(str.c_str()); }

    /// true if 'path' is an existing readible file
    inline bool is_file(std::string const& str) { return is_file(str.c_str()); }

    /// change current working directory
    inline int  change_dir(std::string const& str) { return change_dir(str.c_str()); }

    /// create new directory
    int make_dir(const char path[]);
    
    /// list content of directory
    std::vector<std::string> list_dir(const char path[]);
    
    /// list files in directory `path` that have a name finishing with 'ext'
    std::vector<std::string> list_dir(const char path[], std::string const& ext);
   
    /// extract the directory part from the given path
    std::string dir_part(const char path[]);

    /// extract the file part from the given path
    std::string file_part(const char path[]);

    /// complete the file name using the given directory
    std::string full_name(std::string const& path, std::string const& file);

}

