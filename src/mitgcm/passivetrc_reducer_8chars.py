import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates diagnostic variables of 8 chars
'''
, formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--inputfile','-i',
                                type = str,
                                required = True,
                                help = '''namelist.passivetrc''')
    parser.add_argument(   '--outfile','-o',
                                type = str,
                                required = True,
                                help = '''namelist.passivetrc''')

    return parser.parse_args()

args = argument()



#namelist="/Users/gbolzon/Documents/workspace/ogstm/src/namelists/namelist.passivetrc"
#section="NATTRC_DIAG"
#var_identifier='dianm'

def file2stringlist(filename):
    '''
    Argument : string - a path name indicating a text file
    Returns  : a list of strings

    A text file is converted in a list of strings, one for each row.

    '''
    LIST=[]
    filein=open(filename)
    for line in filein:
        LIST.append(line[:-1])
    filein.close()
    return LIST
def dumpfile(string_list, filename):
    fid = open(filename,'wt')
    for line in string_list:
        fid.write(line+"\n")
    fid.close()

def reducer8(String):
    RULES={"PPY_ii":"_", "BFM1D_exR2ac_ii":"exR2ac", "MEZ_ii":"_",  "MIZ_ii":"_", "OMT_ii":"_", "PBA_ii":"_" }
    for r in RULES.keys():
        String = String.replace(r, RULES[r])
    return String[:8]

F=file2stringlist(args.inputfile)


def stringlist_modify(LINES, section, var_identifier):
    parse=False
    for iline, line in enumerate(LINES):
        if line.strip()==section:
            parse=True
        if parse:
            if var_identifier in line:
                prefix,var, _= line.rsplit("\"")
                new_line = "%s\"%s\"" %(prefix, reducer8(var))
                LINES[iline] = new_line
            if line=="/":
                break
    return LINES

F = stringlist_modify(F, section="&NATTRC_DIAG", var_identifier='dianm')
F = stringlist_modify(F, section="&NATTRC_DIAG_2D", var_identifier='dianm')

dumpfile(F, args.outfile)
       

