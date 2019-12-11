#!/usr/bin/env python3
#
# read_config.py is a simple parser for Cytosim configuration files
#
# F. Nedelec, 13-15.11.2013 - 31.07.2014, 3.09.2015, 12.08.2018, 11.12.2019

"""
    Read cytosim configuration files

Syntax:
   
    read_config.py FILE
    read_config.py - FILE
    
Description:

    Reads a config file and prints a formatted copy to standard output.
    Two formats are possible: vertical (default) or horizontal (second syntax)

F. Nedelec 11.2013--20.02.2017
"""

import sys, os, io

read_comments = False
column_width = 16
err = sys.stderr
out = sys.stdout


def format_values(val):
    """
    Convert values to a cannonical ASCII representation
    """
    if isinstance(val, str):
        return val
    if isinstance(val, dict):
        r = "( "
        for k in sorted(val):
            r += str(k) + '=' + format_values(val[k]) + "; "
        return r + ")"
    try:
        r = format_values(val[0])
        e = len(val) - 1
        if isinstance(val[e], str) and val[e][0] == '%':
            e = e - 1
        for k in val[1:e+1]:
            r += ", " + format_values(k)
        return r
    except:
        return repr(val)


class Instruction:
    words = []
    vals = {}
    cnt  = 1
    
    def __init__(self, key):
        self.words = [ key ]
    
    def value(self, key, index=0):
        if key in self.vals:
            val = self.vals[key]
            if isinstance(val, list):
                return val[index]
            else:
                return val
        return ""
    
    def values(self):
        return self.vals;
    
    def __repr__(self):
        return self.format_horizontal()
    
    def format_header(self):
        res = self.words[0] + ' '
        if self.cnt != 1:
            res += repr(self.cnt) + ' '
        if len(self.words) > 1:
            res += self.words[1] + ' '
        if len(self.words) > 2:
            res += self.words[2] + ' '
        return res
    
    def format_horizontal(self):
        res = self.format_header() + "{ "
        for p in sorted(self.vals):
            res += str(p) + '=' + format_values(self.vals[p]) + "; "
        res += "}"
        return res
    
    def format_vertical(self):
        res = self.format_header() + "\n{\n"
        for p in sorted(self.vals):
            lin = "   " + p.ljust(column_width) + " = " + format_values(self.vals[p])
            #add line-comment if present:
            try:
                e = self.vals[p][-1]
                if e[0] == '%':
                    lin = (lin+';').ljust(50) + e
                else:
                    lin += ';'
            except:
                lin += ';'
            res += lin + '\n'
        res += "}\n"
        return res
    
    def format(self, mode):
        if mode == 1:
            return self.format_horizontal()
        else:
            return self.format_vertical()



def get_char(fid):
    return fid.read(1)


def get_hexadecimal(fid, s):
    res = s
    pos = fid.tell()
    c = get_char(fid)
    while c and c in "0123456789ABCDEFabcdef":
        res += c
        pos = fid.tell()
        c = get_char(fid)
    fid.seek(pos, 0)
    return res


def get_number(fid, s):
    """
    Read a number with decimal point and optional exponent
    """
    res = s
    pos = fid.tell()
    c = get_char(fid)
    if c == 'x':
        return get_hexadecimal(fid, s+c)
    if c == '-' or c == '+':
        res += c
        pos = fid.tell()
        c = get_char(fid)
    no_point = 1
    while c.isdigit() or ( c=='.' and no_point ):
        if c=='.':
            no_point = 0
        res += c
        pos = fid.tell()
        c = get_char(fid)
    if c == 'e' or c == 'E':
        return res + get_number(fid, 'e')
    if c == 'P':
        return res
    fid.seek(pos, 0)
    return res


def delimiter(c):
    if c == '(': return ')'
    if c == '{': return '}'
    if c == '[': return ']'
    if c == '"': return '"'
    return 0


def get_until(fid, e):
    res = ''
    c = get_char(fid)
    while c and c != e:
        res += c
        c = get_char(fid)
    return res


def get_block(fid, s, e):
    """
    Return a block including the enclosing delimiters
    """
    res = s
    c = get_char(fid)
    while c:
        if c == e:
            res += c
            return res
        if delimiter(c):
            res += get_block(fid, c, delimiter(c))
        else:
            res += c
        c = get_char(fid)
    err.write("  Error: unmatched delimiter "+s+"\n")
    sys.exit(1)
    return res


def valid_token_char(c):
    return c.isalnum() or c == ':' or c == '_' or c == '.'


def get_token(fid):
    """
    Extract the next token from the file
    """
    c = get_char(fid)
    while c.isspace() and c != '\n':
        c = get_char(fid)
    if c == '%':
        return '%' + fid.readline()
    if delimiter(c):
        return get_block(fid, c, delimiter(c))
    if c.isdigit() or c == '-' or c == '+':
        return get_number(fid, c)
    res = c
    if valid_token_char(c):
        while c:
            pos = fid.tell()
            c = get_char(fid)
            if not valid_token_char(c):
                fid.seek(pos, 0)
                break
            res += c
    #print(" TOKEN |"+res+"|")
    return res


def uncode(arg):
    try:
        if isinstance(arg, unicode):
            return str(arg.decode('utf-8'))
    except:
        pass
    return arg


def convert(arg):
    """
    Convert a string to a int or float if possible
    """
    try:
        return int(arg)
    except:
        pass
    try:
        return float(arg)
    except:
        pass
    #remove spaces:
    s = arg.strip()
    #remove outer delimiters:
    if len(s) > 3 and delimiter(s[0]) == s[-1]:
        s = s[1:-1].strip()
    return s


def simplify(arg):
    """
    Simplify the values of imbricated sets and dictionaries
    """
    #print('SIMPLIFY   ', type(arg), arg)
    arg = uncode(arg)
    if isinstance(arg, str):
        return convert(arg)
    if isinstance(arg, dict):
        s = {}
        for k in arg:
            uk = uncode(k)
            s[uk] = simplify(arg[k])
        #print('SIMPLIFY=> ', s)
        return s
    try:
        if len(arg) == 1:
            return simplify(arg[0])
        s = []
        for a in arg:
            s.append(simplify(a))
        #print('SIMPLIFY=> ', s)
        return s
    except:
        return arg


def file_object(arg):
    #print("new_file_object", type(arg))
    try:
        return io.StringIO(arg)
    except:
        return io.StringIO(unicode(arg))


def read_list(fid):
    """
    Process a list of key=values
    """
    dic = {}
    key = ''
    val = ''
    kkk = ''
    tok = get_token(fid)
    while tok:
        #print('key %10s token |%s|' %(key, tok))
        if tok == '=':
            key = kkk
            dic[key] = []
            val = ''
        elif tok == '\n' or tok == ';' or tok[0] == '%':
            if key:
                dic[key].append(val)
                val = ''
                key = ''
            if tok[0] == '%':
                global read_comments
                if read_comments and kkk and len(tok) > 1:
                    dic[kkk].append(tok)
            if not tok == ';':
                kkk = ''
        elif tok == ',':
            if key:
                dic[key].append(val)
                val = ''
        elif delimiter(tok[0]):
            val = read_list(file_object(tok[1:-1]))
            if not val:
                val = tok
        else:
            if key:
                if val:
                    val += ' ' + tok
                else:
                    val = tok
            else:
                kkk = tok
        tok = get_token(fid)
    if key and val:
        dic[key].append(val)
    return dic


def parse_config(fid):
    """
    return the list resulting from parsing the specified file
    """
    known_keywords = [ 'set', 'new', 'run', 'call', 'delete', 'change', 'cut', 'mark' ];
    cur = []
    pile = []
    while fid:
        tok = get_token(fid)
        #print(">>>>> `",tok.strip(), "'")
        if not tok:
            break
        if tok[0] == '%':
            pass
        elif tok == '\n':
            pass
        elif tok == 'repeat':
            if cur:
                pile.append(cur)
                cur = []
            cnt = get_token(fid)
            blok = get_token(fid)
            sfid = io.StringIO(blok[1:-1])
            pile.extend(parse_config(sfid))
        elif tok in known_keywords:
            if cur:
                pile.append(cur)
            cur = Instruction(tok)
        elif len(cur.words)==1 and tok.isdigit():
            cur.cnt = int(tok)
        elif tok[0] == '{' or tok[0] == '(':
            pam = read_list(file_object(tok[1:-1]))
            #print("VALUES  : ", pam)
            cur.vals = simplify(pam)
            #print("SIMPLIFIED: ", cur.vals)
            pile.append(cur)
            cur = []
        elif tok[0].isalpha():
            cur.words.append(tok)
        else:
            print(">>>>>ignored `",tok, "'")
    return pile


def parse(arg):
    """
    Returns a dictionnary obtained by parsing the file specified by name `arg`
    """
    if isinstance(arg, str):
        with open(arg, 'r') as f:
            pile = parse_config(f)
            f.close()
    else:
        pile = parse_config(arg)
    return pile

#------------------------------------------------------------------------

def format_config(pile, mode=0, prefix=''):
    """
    Returns a compact representation of the hierachical dictionnary
    """
    res = ''
    for i in pile:
        res += prefix + i.format(mode) + "\n"
    return res


def get_command(pile, keys):
    """
    Return all commands corresponding to = [ command, class, name ]
    The keys can be '*' to specify a wildcard matching
    """
    if len(keys) != 3:
        return "Error: get_command(arg, keys) expects `keys` to be of size 3"
    res = []
    for p in pile:
        if keys[0]=='*' or p.words[0]==keys[0]:
            if keys[1]=='*' or p.words[1]==keys[1]:
                if keys[2]=='*' or p.words[2]==keys[2]:
                    res.append(p)
    if len(res) == 1:
        return res[0]
    return res


def get_value(arg, keys):
    """
    Return the value speficied by keys = [ command, class, name, parameter ]
    The first 3 keys can be '*' for wildcard matching
    """
    if len(keys) != 4:
        return "Error: get_value(arg, keys) expects `keys` to be of size 4"
    if isinstance(arg, str):
        pile = parse(arg)
    else:
        pile = arg
    for p in pile:
        if keys[0]=='*' or p.words[0]==keys[0]:
            if keys[1]=='*' or p.words[1]==keys[1]:
                if keys[2]=='*' or p.words[2]==keys[2]:
                    try:
                        return p.vals[keys[3]]
                    except KeyError:
                        if keys[3]=='*':
                            return p.vals
    return "Unspecified"


#------------------------------------------------------------------------

def main(args):
    global read_comments;
    horizontal = False
    files = []
    
    for arg in args:
        if os.path.isfile(arg):
            files.append(arg)
        elif arg=="-" or arg=="-h":
            horizontal = True
        elif arg=="+":
            read_comments = True
        else:
            err.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    
    if not files:
        err.write("  Error: you must specify a file to read\n")
        sys.exit()

    for f in files:
        out.write("\n")
        pile = parse(f)
        out.write(format_config(pile, horizontal))
        #print(f.rjust(40) + " " + value(pile, ['set', 'space', '*', 'geometry']) )


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])
