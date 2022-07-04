import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Rewrites BFM_var_list.h and BFM1D_Output_Ecology.F90 using "dump" attribute in BFMtab.xml
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--inputdir','-i',
                                type = str,
                                required = True,
                                default = "CODE/bfm/src/ogstm/ .",
                                help = 'directory whith bfmv5 generated files')
    parser.add_argument(   '--outdir','-o',
                                type = str,
                                required = True,
                                default = "/some/path",
                                help = 'Directory to put reduced BFM_var_list.h and BFM1D_Output_Ecology.F90')
    parser.add_argument(   '--xmlfile','-f',
                                type = str,
                                required = True,
                                default = "CODE/ogstm/bfmv5/BFMtab.xml",
                                help = 'Path of the BFMtab with highfreq and dump info')
    return parser.parse_args()

args = argument()
from xml.dom import minidom
import os,sys
import numpy as np
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
def addsep(string):
    if string[-1] != os.sep:
        return string + os.sep
    else:
        return  string
def lines_to_vars_table(string_list,nvars):
    VARS=np.zeros((nvars,),dtype=[("var","S50"), ("ind",np.int)])
    ivar=0
    for line in string_list:
        if line.find("::")==-1: continue
        prefix, Vars=line.rsplit("::")
        raw_vars=Vars.rsplit(",")
        for raw_var in raw_vars:
            left_side, right_side=raw_var.rsplit("=")
            left_side=left_side.replace(" ","") # deblank
            VARS[ivar]["var"]=left_side[2:]
            VARS[ivar]["ind"]=right_side
            ivar=ivar+1
    return VARS
def varstable_to_lines(VARS_TABLE):
    LINES=[]
    for ivar, var in enumerate(VARS_TABLE['var']):
        s="     integer,parameter ::  pp%s=%d ! ogstm reduced list index" %(var, ivar+1)
        LINES.append(s)
    return LINES
def vars_table_to_lines_Output_Ecology(string_list, VARS_TABLE,ndims=3):
    if ndims==3:
        string="local_BFM1D_dia("
        suffix=",:"
    if ndims==2:
        string="local_BFM0D_dia2d("
        suffix=")"
    LINES_RED=[]
    for ivar, var in enumerate(VARS_TABLE['var']):
        orig_index= VARS_TABLE[ivar]['ind']
        str_to_find="pp" + var
        pos=str_to_find.find("_iiP")
        if pos >-1 :
            str_to_find = str_to_find[:pos] + "("  + str_to_find[pos+1:] + ")"

        for line in string_list:
            if line.find(str_to_find)>-1:
                left_side,right_side=line.rsplit(string)
                pos_sep = right_side.find(suffix)
                new_line = "%s%d%s" %(string, ivar+1, right_side[pos_sep:])
                LINES_RED.append(new_line)
    return LINES_RED



def xml_read(filename):
    xmldoc = minidom.parse(filename)
    NODES=xmldoc.getElementsByTagName("Diagnostics")[0].getElementsByTagName("var")
    XML_MODELVARS={}
    for node in NODES:
        var = str(node.getAttribute("name"))
        high= node.getAttribute("dump")=='true'

        if not XML_MODELVARS.has_key(var): # check about BFMtab formatting
            XML_MODELVARS[var] = high
        else:
            print var + " is already defined in xml"
            sys.exit(1)
    NODES=xmldoc.getElementsByTagName("Diagnostics_2D")[0].getElementsByTagName("var")
    for node in NODES:
        var = str(node.getAttribute("name"))
        high= node.getAttribute("dump")=='true'

        if not XML_MODELVARS.has_key(var): # check about BFMtab formatting
            XML_MODELVARS[var] = high
        else:
            print var + " is already defined in xml"
            sys.exit(1)
    return XML_MODELVARS

def dumpfile(string_list, filename):
    fid = open(filename,'wt')
    for line in string_list:
        fid.write(line+"\n")
    fid.close()


INPUTDIR=addsep(args.inputdir)
OUTDIR  =addsep(args.outdir)
LINES=file2stringlist(INPUTDIR + "BFM_var_list.h")
XML_MODELVARS = xml_read(args.xmlfile)
nLINES=len(LINES)

for iline, line in enumerate(LINES):
    if line.find("INTEGER, parameter :: jptra = ")>-1:
        i_jptra_line=iline

    if line.find("INTEGER, parameter :: jptra_var")>-1:
        left_side,right_side = line.rsplit("INTEGER, parameter :: jptra_var =")
        jptra_var = int(right_side.replace(" ",""))

    if line.find("INTEGER, parameter :: jptra_flux")>-1:
        left_side,right_side = line.rsplit("INTEGER, parameter :: jptra_flux =")
        jptra_flux = int(right_side.replace(" ",""))

    if line.find("INTEGER, parameter :: jptra_dia_2d")>-1:
        i_jptra_line_2d=iline
        left_side,right_side = line.rsplit("INTEGER, parameter :: jptra_dia_2d =")
        jptra_dia_2d = int(right_side.replace(" ",""))

    if line.find("diagnostic indexes")>-1:
        line_start=iline+1
    if line.find("flux indexes")>-1:
        line_end=iline+2 # in order to read following line
    if line.find("variables 2d")>-1:
        line_var2d=iline

jptra_dia = jptra_var + jptra_flux
DIAGNOSTIC_LINES_3D  =LINES[line_start:line_end]
DIAGNOSTIC_LINES_2D  =LINES[line_var2d:]

VARS3D = lines_to_vars_table(DIAGNOSTIC_LINES_3D, jptra_dia)
VARS2D = lines_to_vars_table(DIAGNOSTIC_LINES_2D, jptra_dia_2d)

good3D=np.zeros((jptra_dia   ,),np.bool)
good2D=np.zeros((jptra_dia_2d,),np.bool)

for ivar in range(jptra_dia):
    varname=VARS3D[ivar]["var"]
    good3D[ivar]=XML_MODELVARS[varname]
for ivar in range(jptra_dia_2d):
    varname=VARS2D[ivar]["var"]
    good2D[ivar]=XML_MODELVARS[varname]

jptra_dia_red_3d=good3D.sum()
jptra_dia_red_2d=good2D.sum()

VARS3D_TO_DUMP=VARS3D[good3D]
VARS2D_TO_DUMP=VARS2D[good2D]

DIAGNOSTIC_LINES_RED_3D= varstable_to_lines(VARS3D_TO_DUMP)
DIAGNOSTIC_LINES_RED_2D= varstable_to_lines(VARS2D_TO_DUMP)

#BFM_var_list.h reconstruction
LINES_RED=[]

LINES_RED.append(LINES[i_jptra_line])
LINES_RED.append("")
LINES_RED.append(("      INTEGER, parameter :: jptra_dia = %d" % (jptra_dia_red_3d)))
LINES_RED.append("")
LINES_RED.append(("      INTEGER, parameter :: jptra_dia_2d = %d" % (jptra_dia_red_2d)))

for i in range(i_jptra_line_2d+1,line_start): LINES_RED.append(LINES[i]) # state
for line in DIAGNOSTIC_LINES_RED_3D: LINES_RED.append(line)              # dia3D
for i in range(line_end,line_var2d+1): LINES_RED.append(LINES[i])           # flux
for line in DIAGNOSTIC_LINES_RED_2D: LINES_RED.append(line)


dumpfile(LINES_RED, OUTDIR + "BFM_var_list.h")


inputfile=INPUTDIR + "BFM1D_Output_Ecology.F90"
LINES=file2stringlist(inputfile)
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("INTEGER, parameter :: jptra_var")>-1:
        i_jptra_var_line=iline
    if line.find("INTEGER, parameter :: jptra_dia_2d")>-1:
        i_jptra_var_line_2d=iline
    if line.find("  local_BFM1D_dia = 0;")> -1:
        line_start=iline
    if line.find("local_BFM0D_dia2d = 0;")> -1:
        line_mid=iline
    if line.find("! ***********************************  First Group Done")> -1: # prende il secondo
        line_end=iline

DIAGNOSTIC_LINES_3D=LINES[line_start+3:line_mid-1]
DIAGNOSTIC_LINES_2D=LINES[line_mid:line_end]
DIAGNOSTIC_LINES_RED_3D = vars_table_to_lines_Output_Ecology(DIAGNOSTIC_LINES_3D, VARS3D_TO_DUMP, ndims=3)
DIAGNOSTIC_LINES_RED_2D = vars_table_to_lines_Output_Ecology(DIAGNOSTIC_LINES_2D, VARS2D_TO_DUMP, ndims=2)

#BFM1D_Output_Ecology.F90 reconstruction
LINES_RED=[]

for i in range(i_jptra_var_line_2d+1):
    if i==i_jptra_var_line:
        LINES_RED.append("      INTEGER, parameter :: jptra_var = %d" % (jptra_dia_red_3d))
    else:
        if i==i_jptra_var_line_2d:
            LINES_RED.append("      INTEGER, parameter :: jptra_dia_2d = %s" %(jptra_dia_red_2d))
        else:
            LINES_RED.append(LINES[i])

for i in range(i_jptra_var_line_2d+1,line_start+2): LINES_RED.append(LINES[i])
for line in DIAGNOSTIC_LINES_RED_3D: LINES_RED.append(line)
for i in range(line_mid,line_mid+2): LINES_RED.append(LINES[i])
for line in DIAGNOSTIC_LINES_RED_2D: LINES_RED.append(line)
for i in range(line_end,nLINES) : LINES_RED.append(LINES[i])
dumpfile(LINES_RED, OUTDIR + "BFM1D_Output_Ecology.F90")

