#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# MODEL  BFM - Biogeochemical Flux Model
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#
# ROUTINE: diff_apply.py
#
# DESCRIPTION
#    Writes code for MITgcm coupler
#
# COPYING
#
#   Copyright (C) 2023 BFM System Team (bfm_st@cmcc.it)
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation.
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#   See the GNU General Public License for more details.
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
import os
import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
   Writes:
   gchem_init_vari.F
   gchem_init_fixed.F
   gchem_fields_load.F
   gchem_readparms.F
   gchem_calc_tendency.F
   longstep_gchem_calc_tendency.F"
   GCHEM.h
   GCHEM_OPTIONS.h
   longstep_thermodynamics.F
   ptracers_reset.F
   GCHEM_FIELDS.H

   by reading from MITgcm code
'''
, formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--inputdir','-i',
                                type = str,
                                required = True,
                                help = '''namelist.passivetrc''')
    parser.add_argument(   '--outdir','-o',
                                type = str,
                                required = True,
                                help = '''path of the output dir''')

    return parser.parse_args()

args = argument()



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



def insert_lines(orig_lines,NEW_LINES,position_line,nLINES,final=False):
    OUTLINES=[]
    for iline in range(position_line):
        OUTLINES.append(orig_lines[iline])

    for line in NEW_LINES: OUTLINES.append(line)

    for iline in range(position_line,nLINES):
        OUTLINES.append(orig_lines[iline])
    if final:
        OUTLINES=[line + "\n" for line in OUTLINES]
    return OUTLINES
def replace_lines(orig_lines, old, new_lines):
    OUTLINES=[]
    found=False
    for line in orig_lines:
        if line.find(old) > -1 :
            found = True
            for dest_line in new_lines:
                OUTLINES.append(dest_line + "\n")
        else:
            OUTLINES.append(line)
    if not found :
        print(old, "Not found")
        return None
    return OUTLINES

INPUTDIR=addsep(args.inputdir)
OUTDIR=addsep(args.outdir)
MITCODE=INPUTDIR + "pkg/gchem/"
MYCODE=OUTDIR


filename="gchem_init_vari.F"
infile=MITCODE + filename
outfile=MYCODE + filename

LINES=file2stringlist(infile)
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("INTERFACE: ==") >-1: position_line=iline
   
NEW_LINES=[
"#ifdef ALLOW_BFMCOUPLER",
"#include \"BFMcoupler_OPTIONS.h\"",
"#endif"]
OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES)

LINES=OUTLINES
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("#endif /* ALLOW_GCHEM */") >-1: position_line=iline
NEW_LINES=[
"#ifdef ALLOW_BFMCOUPLER",
"      IF ( useBFMcoupler) THEN",
"         CALL BFMcoupler_INI_FORCING(myThid)",
"      ENDIF",
"#endif"]
OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES,final=True)
fid=open(outfile,'w')
fid.writelines(OUTLINES)
fid.close()




filename="gchem_init_fixed.F"
infile=MITCODE + filename
outfile=MYCODE + filename
LINES=file2stringlist(infile)
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("#ifdef ALLOW_DIAGNOSTICS") >-1: position_line=iline
NEW_LINES=[
"#ifdef ALLOW_BFMCOUPLER",
"         call BFMcoupler_INIT_FIXED(myThid)",
"#endif"]
OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES,final=True)
fid=open(outfile,'w')
fid.writelines(OUTLINES)
fid.close()


filename="gchem_fields_load.F"
infile=MITCODE + filename
outfile=MYCODE + filename
LINES=file2stringlist(infile)
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("#endif /* ALLOW_GCHEM */") >-1: position_line=iline
NEW_LINES=[
"#ifdef ALLOW_BFMCOUPLER",
"      IF ( useBFMcoupler ) THEN",
"       CALl BFMcoupler_FIELDS_LOAD(myIter,myTime,myThid)",
"      ENDIF",
"#endif /* ALLOW_BFMCOUPLER */"]
OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES,final=True)
fid=open(outfile,'w')
fid.writelines(OUTLINES)
fid.close()


filename="gchem_readparms.F"
infile=MITCODE + filename
outfile=MYCODE + filename

LINES=file2stringlist(infile)
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("NAMELIST /GCHEM_PARM01/") >-1: position_line=iline+1
NEW_LINES=[
"#ifdef ALLOW_BFMCOUPLER",   
"c                 useBFMcoupler must be read in namelist",
"     &            useBFMcoupler,",
"#endif"]
OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES)

LINES=OUTLINES
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("C Set defaults values for parameters in GCHEM.h") >-1: position_line=iline+1
NEW_LINES=[
"#ifdef ALLOW_BFMCOUPLER",
"       useBFMcoupler = .FALSE.",
"#endif"]
OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES)

LINES=OUTLINES
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("#endif /* ALLOW_GCHEM */")>-1: position_line=iline
NEW_LINES=[
"#ifdef ALLOW_BFMCOUPLER",
"      IF ( useBFMcoupler ) THEN",
"        CALL BFMcoupler_READPARMS(myThid)",
"      ENDIF",
"#endif /* ALLOW_BFMCOUPLER */"]     
OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES,final=True)

fid=open(outfile,'w')
fid.writelines(OUTLINES)
fid.close()


filename="gchem_calc_tendency.F"
infile=MITCODE + filename
outfile=MYCODE + filename

LINES=file2stringlist(infile)
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("# endif /* GCHEM_SEPARATE_FORCING */")>-1: position_line=iline
NEW_LINES=[
"C------------------------",
"C BFM coupler           |",
"C------------------------",
"#ifdef ALLOW_BFMCOUPLER",
"      IF ( useBFMcoupler ) THEN",
"C$taf loop = parallel",
"        DO bj=myByLo(myThid),myByHi(myThid)",
"        DO bi=myBxLo(myThid),myBxHi(myThid)",
"              ",
"        jMin=1",
"        jMax=sNy",
"        iMin=1",
"        iMax=sNx",
"C BFMcoupler operates on bi,bj part only, but needs to get full arrays i",
"C because of last index (iPtr)",
"         CALL BFMcoupler_calc_tendency(bi,bj,imin,imax,jmin,jmax,",
"     &                            myIter,myTime,myThid)",
"        ENDDO",
"        ENDDO",
"       ENDIF",
"#endif /* ALLOW_BFMCOUPLER */",
"#endif /* ALLOW_LONGSTEP */"]


OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES,final=False)
LINES=OUTLINES
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("# ifndef GCHEM_SEPARATE_FORCING")>-1: position_line=iline
NEW_LINES=["# ifndef ALLOW_LONGSTEP"]
nLINES = len(LINES)
OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES,final=True)
fid=open(outfile,'w')
fid.writelines(OUTLINES)
fid.close()


outfile=MYCODE + 'longstep_gchem_calc_tendency.F'
LINES=OUTLINES # we'll apply changes from gchem_calc_tendency.F

l0="C $Header: /u/gcmpack/MITgcm/pkg/gchem/gchem_calc_tendency.F,v 1.5 2013/06/10 02:52:57 jmc Exp $"
l1="C $Header: longstep_gchem_calc_tendency.F,v 1.5 2015/04/01 02:52:57 GPC Exp $"
OUTLINES = replace_lines(LINES,l0,[l1])

l0="C !ROUTINE: GCHEM_CALC_TENDENCY"
l1="C !ROUTINE: LONGSTEP_GCHEM_CALC_TENDENCY"
OUTLINES = replace_lines(OUTLINES,l0,[l1] )


l0= "      SUBROUTINE GCHEM_CALC_TENDENCY("
l1 = ["C version of GCHEM_CALC_TENDENCY called by LONGSTEP package",
 "      SUBROUTINE LONGSTEP_GCHEM_CALC_TENDENCY("]
OUTLINES = replace_lines(OUTLINES,l0, l1 )

l0 = "# ifndef ALLOW_LONGSTEP"
l1 = "# ifdef ALLOW_LONGSTEP"
OUTLINES = replace_lines(OUTLINES,l0,[l1] )

fid=open(outfile,'w')
fid.writelines(OUTLINES)
fid.close()





filename="GCHEM.h"
infile=MITCODE + filename
outfile=MYCODE + filename
LINES=file2stringlist(infile)
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("C     useDARWIN :: flag to turn on/off darwin pkg")>-1: position_line=iline+1
NEW_LINES=[
"#ifdef ALLOW_BFMCOUPLER",
"C     useBFMcoupler :: flag to turn on/off BFMcoupler pkg",
"#endif"]
OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES)
LINES=OUTLINES
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("     &              useDARWIN")>-1: position_line=iline
NEW_LINES=[
"#ifdef ALLOW_BFMCOUPLER",
"     &              ,useBFMcoupler",
"      LOGICAL useBFMcoupler",
"#endif"]
OUTLINES=insert_lines(LINES, NEW_LINES, position_line+1,nLINES,final=True)
fid=open(outfile,'w')
fid.writelines(OUTLINES)
fid.close()


filename="GCHEM_OPTIONS.h"
infile=MITCODE + filename
outfile=MYCODE + filename
LINES=file2stringlist(infile)
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("#endif /* ALLOW_GCHEM */")>-1: position_line=iline
NEW_LINES=[
"#undef GCHEM_SEPARATE_FORCING",
"c undefining gchem_separate_forcing actives BFMcoupler_calc_tendency and add_tendency",
"c  #define GCHEM_SEPARATE_FORCING"]    
OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES,final=True)
fid=open(outfile,'w')
fid.writelines(OUTLINES)
fid.close()    
    

##### applying differences in longstep Pkg
MITCODE=INPUTDIR + "pkg/longstep/"
filename="longstep_thermodynamics.F"
infile=MITCODE + filename
outfile=MYCODE + filename
LINES=file2stringlist(infile)
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("      IF ( LS_doTimeStep ) THEN")>-1: position_line=iline+1
NEW_LINES=[
"cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc",
"c CGP 2015/04/03 adding call to gchem_calc_tendency",
"      CALL LONGSTEP_GCHEM_CALC_TENDENCY( myTime, myIter, myThid )",
"cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"]
OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES,final=True)
fid=open(outfile,'w')
fid.writelines(OUTLINES)
fid.close()    
    

##### applying differences in ptracer Pkg
MITCODE=INPUTDIR + "pkg/ptracers/"
filename="ptracers_reset.F"
infile=MITCODE + filename
outfile=MYCODE + filename
LINES=file2stringlist(infile)
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("C     !INPUT PARAMETERS:")>-1: position_line=iline
NEW_LINES=[
"#ifdef ALLOW_GCHEM",
"#include \"GCHEM.h\"",
"#endif"]  
OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES)
LINES=OUTLINES
nLINES=len(LINES)
for iline, line in enumerate(LINES):
    if line.find("     DO iTracer = 1, PTRACERS_num")>-1:position_line=iline+1 
NEW_LINES=[
"c check for negative values of pTracer variables and set them to 1._d-10",
"#ifdef ALLOW_GCHEM",
"      IF ( useGCHEM ) THEN",
"#ifdef ALLOW_BFMCOUPLER",
"         IF (useBFMcoupler) THEN",
"           DO bj = myByLo(myThid), myByHi(myThid)",
"             DO bi = myBxLo(myThid), myBxHi(myThid)",
"               DO k=1,Nr",
"                 DO j=1-OLy,sNy+OLy",
"                   DO i=1-OLx,sNx+OLx",
"                     if (pTracer(i,j,k,bi,bj,iTracer).lt.0.0)THEN",
"                         pTracer(i,j,k,bi,bj,iTracer)=1. _d -10",
"                     ENDIF",
"                   ENDDO",
"                 ENDDO",
"               ENDDO",
"             ENDDO",
"           ENDDO",
"         ENDIF",
"#endif /* BFMCOUPLER */",
"      ENDIF",
"#endif /* ALLOW_GCHEM */ "]
OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES,final=True)
fid=open(outfile,'w')
fid.writelines(OUTLINES)
fid.close()    



filename="GCHEM_FIELDS.h"
MITCODE=INPUTDIR + "pkg/gchem/"
infile=MITCODE + filename
outfile=MYCODE + filename
LINES=file2stringlist(infile)
nLINES = len(LINES)

for iline, line in enumerate(LINES):
    if line.find("#endif /* GCHEM_SEPARATE_FORCING */")>-1: position_line=iline
NEW_LINES=["#ifdef GCHEM_SEPARATE_FORCING",
"      _RL gchemTendency(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,",
"     &                  PTRACERS_num)",
"      COMMON /GCHEM_FIELDS/",
"     &     gchemTendency",
"#endif /* when define GCHEM_SEPARATE_FORCING */"]
OUTLINES=insert_lines(LINES, NEW_LINES, position_line,nLINES,final=True)
fid=open(outfile,'w')
fid.writelines(OUTLINES)
fid.close()

















