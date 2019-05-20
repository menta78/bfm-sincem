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

INPUTDIR=addsep(args.inputdir)
OUTDIR  =addsep(args.outdir)
LINES=file2stringlist(INPUTDIR + "BFM_var_list.h")
nLINES=len(LINES)

for iline, line in enumerate(LINES):
    if line.find("INTEGER, parameter :: jptra_var")>-1:
        i_jptra_line=iline
        left_side,right_side = line.rsplit("INTEGER, parameter :: jptra_var =")
        jptra_dia = int(right_side.replace(" ",""))
        
    if line.find("diagnostic indexes")>-1:
        line_start=iline+1
    if line.find("flux indexes")>-1:
        line_end=iline-1

DIAGNOSTIC_LINES=LINES[line_start:line_end]

VARS=np.zeros((jptra_dia,),dtype=[("var","S50"), ("ind",np.int)])
good=np.zeros((jptra_dia,),np.bool)
ivar=0
for line in DIAGNOSTIC_LINES: 
    #print line
    if line.find("::")==-1: continue
    prefix, Vars=line.rsplit("::")
    raw_vars=Vars.rsplit(",")
    for raw_var in raw_vars:
        left_side, right_side=raw_var.rsplit("=")
        left_side=left_side.replace(" ","") # deblank
        VARS[ivar]["var"]=left_side[2:]
        VARS[ivar]["ind"]=right_side
        ivar=ivar+1

xmldoc = minidom.parse(INPUTDIR + 'BFMtab.xml')
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

for ivar in range(jptra_dia):
    varname=VARS[ivar]["var"]
    good[ivar]=XML_MODELVARS[varname]

jptra_dia_red=good.sum()
VARS_TO_DUMP=VARS[good]
DIAGNOSTIC_LINES_RED=[]
for ivar, var in enumerate(VARS_TO_DUMP['var']):
    s="integer,parameter ::  pp%s=%d" %(var, ivar+1)
    DIAGNOSTIC_LINES_RED.append(s)

#BFM_var_list.h reconstruction
LINES_RED=[]

for i in range(i_jptra_line): LINES_RED.append(LINES[i])
LINES_RED.append("      INTEGER, parameter :: jptra_var = %d" % (jptra_dia_red))
for i in range(i_jptra_line+1,line_start): LINES_RED.append(LINES[i])
for line in DIAGNOSTIC_LINES_RED: LINES_RED.append(line)
for i in range(line_end,nLINES): LINES_RED.append(LINES[i])

outfile=OUTDIR + "BFM_var_list.h"
fid = open(outfile,'wt')
for line in LINES_RED:
    fid.write(line+"\n")
fid.close()


inputfile=INPUTDIR + "BFM1D_Output_Ecology.F90"
outfile  =OUTDIR + "BFM1D_Output_Ecology.F90"
LINES=file2stringlist(inputfile)
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("INTEGER, parameter :: jptra_var")>-1:
        i_jptra_line=iline
    if line.find("  local_BFM1D_dia = 0;")> -1:
        line_start=iline+2
    if line.find("local_BFM0D_dia2d = 0;")> -1:
        line_end=iline-1

DIAGNOSTIC_LINES=LINES[line_start+1:line_end]
DIAGNOSTIC_LINES_RED=[]
for ivar, var in enumerate(VARS_TO_DUMP['var']):
    orig_index= VARS_TO_DUMP[ivar]['ind']
    str_to_find="local_BFM1D_dia(%d,:)" % orig_index
    for line in DIAGNOSTIC_LINES:
        if line.find(str_to_find)>-1:
            left_side,right_side=line.rsplit("local_BFM1D_dia(")
            pos_sep = right_side.find(",:")
            new_line = "local_BFM1D_dia(%d%s" %(ivar, right_side[pos_sep:])
            DIAGNOSTIC_LINES_RED.append(new_line)

#BFM1D_Output_Ecology.F90 reconstruction
LINES_RED=[]

for i in range(i_jptra_line): LINES_RED.append(LINES[i])
LINES_RED.append("      INTEGER, parameter :: jptra_var = %d" % (jptra_dia_red))
for i in range(i_jptra_line+1,line_start): LINES_RED.append(LINES[i])
for line in DIAGNOSTIC_LINES_RED: LINES_RED.append(line)
for i in range(line_end,nLINES): LINES_RED.append(LINES[i])

fid = open(outfile,'wt')
for line in LINES_RED:
    fid.write(line+"\n")
fid.close()
