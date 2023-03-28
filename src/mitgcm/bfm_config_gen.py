import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
BFMcoupler configuration script

The script parses the 'namelist.passivetrc' bfm parameter file
in order to generate the needed BFMcoupler header files accordingly.

Main class is bfm_vars, which contains the main vars list
and provides all the methods that need to be called
for the parameter files generation
'''
, formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--inputfile','-i',
                                type = str,
                                required = True,
                                help = '''namelist.passivetrc''')
    parser.add_argument(   '--type','-t',
                                type = str,
                                required = True,
                                choices = ['code','namelist'],
                                help = '''Flag to generate code or namelists''')
    parser.add_argument(   '--outdir','-o',
                                type = str,
                                required = True,
                                help = '''path of the output dir''')

    return parser.parse_args()

args = argument()

import re, os

def addsep(string):
    if string[-1] != os.sep:
        return string + os.sep
    else:
        return  string

OUTDIR=addsep(args.outdir)


def parse_namelist(namelist, section, var_identifier, units_identifier):

    var = list()
    units = list()
    parse = False

    with open(namelist, 'r') as nml:
        for line in nml:
            if ('&' + section) == line.strip():
                parse = True
            elif line.strip() == '/':
                parse = False
            elif parse:
                if var_identifier in line:
                    val = re.findall(r'\"([^"]*)\"', line)[0]
                    var.append(val)
                elif units_identifier in line:
                    val = re.findall(r'\"([^"]*)\"', line)[0]
                    units.append(val)
    n=len(var)
    return [(var[k], units[k]) for k in range(n)]


class BFM_vars(object):



    def __init__(self, namelist,
                 var_section='NATTRC', var_identifier='ctrcnm', var_units_identifier='ctrcun',
                 diag_var_section='NATTRC_DIAG', diag_var_identifier='dianm', diag_var_units_identifier='diaun',
                 diag_var_section_2d='NATTRC_DIAG_2D', diag_var_identifier_2d='dianm', diag_var_units_identifier_2d='diaun'):

        self._vars = parse_namelist(namelist, var_section, var_identifier, var_units_identifier)
        self._vars_properties = dict() # additional properties for all tracers, initialized to null
        diag_vars = parse_namelist(namelist, diag_var_section, diag_var_identifier, diag_var_units_identifier)
        diag_vars_2d = parse_namelist(namelist, diag_var_section_2d, diag_var_identifier_2d, diag_var_units_identifier_2d)
        self._diag_vars = diag_vars + diag_vars_2d # diag_vars and diag_vars_2d in one single list
        self._diag_vars_MIT = list() # MITgcm-BFM-only diagnostic variables, initialized to null



    def add_var(self, var, units):
        self._vars.append((var, units))



    def add_diag_var(self, var, units):
        self._diag_vars.append((var, units))



    def add_diag_var_MIT(self, var, units, idx): # MIT-only diagnostic variables have an additional index
        self._diag_vars_MIT.append((idx, (var, units)))



    def set_vars_property(self, name, value):
        self._vars_properties[name] = value



    def flush(self):

        print('Vars:\n')
        print self._vars
        print('\nDiagnostic variables:\n')
        print self._diag_vars
        print('\nDiagnostic variables - MITgcm-BFM-only:\n')
        print self._diag_vars_MIT



    def write_diagnostic_init(self, output_file='BFMcoupler_diagnostics_init.F'):

        # First part of the code
        header = '''C $Header:/MITgcm/pkg/BFMcoupler/BFMcoupler_diagnostics_init.F,v 1.01
C 2014/04/10

#include "GCHEM_OPTIONS.h"

C !INTERFACE: ==========================================================
      SUBROUTINE BFMcoupler_DIAGNOSTICS_INIT(myThid )

C !DESCRIPTION:
C define diagnostics for BFMcoupler package
C experiment

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS

C     !LOCAL VARIABLES:
      INTEGER        diagNum
      CHARACTER*8    diagName
      CHARACTER*16   diagCode
      CHARACTER*16   diagUnits
      CHARACTER*(80) diagTitle

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Define diagnostics Names :
'''

        # diag_code(1)='S'; % Scalar
        # diag_code(2)='M'; % C-Grid mass point
        # diag_code(3)='R'; % same but cumulate product by hFac & level thickness
        # diag_code(4)='P'; % positive definite diagnostic
        # diag_code(5-8) not used
        # diag_code(9)='M'; % model-level middle
        # diag_code(10)='R'; % levels = Nr#
        diag_code = 'SMRP    MR'

        # Last part of the code
        footer = '''C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END'''

        # write to file
        with open(output_file, 'w') as ofile:

            ofile.write(header)

            # loop on diagnostic variables
            for var in self._diag_vars:
                ofile.write('      diagName  = \'' + var[0] + '\'\n')
                ofile.write('      diagTitle = \'' + var[0] + '\'\n') # long name set equal to name
                ofile.write('      diagUnits = \'' + var[1] + '\'\n')
                ofile.write('      diagCode  = \'' + diag_code + '\'\n')
                ofile.write('      CALL DIAGNOSTICS_ADDTOLIST( diagNum,' + '\n')
                ofile.write('     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )' + '\n')
                ofile.write('C ' + '\n')

            # loop on MITgcm-BFM-only diagnostic variables
            for var in self._diag_vars_MIT:
                ofile.write('      diagName  = \'' + var[1][0] + '\'\n')
                ofile.write('      diagTitle = \'' + var[1][0] + '\'\n') # long name set equal to name
                ofile.write('      diagUnits = \'' + var[1][1] + '\'\n')
                ofile.write('      diagCode  = \'' + diag_code + '\'\n')
                ofile.write('      CALL DIAGNOSTICS_ADDTOLIST( diagNum,' + '\n')
                ofile.write('     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )' + '\n')
                ofile.write('C ' + '\n')

            ofile.write(footer)



    def write_local(self, output_file='BFMcoupler_VARDIAGlocal.h'):

        # First part of the code
        header = '''C $Header:/MITgcm/pkg/BFMcoupler/BFMcoupler_VARDIAGlocal.h,v 1.01
C 2014/04/10
C local variables of diagnostic for BFMcoupler_calc_tendency.f
C 
'''

        # write to file
        with open(output_file, 'w') as ofile:

            ofile.write(header)

            # loop on diagnostic variables
            for var in self._diag_vars:
                ofile.write('      _RL   dia' + var[0] + '(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)\n')

            # loop on MITgcm_BFM-only diagnostic variables
            for var in self._diag_vars_MIT:
                ofile.write('      _RL   dia' + var[1][0] + '(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)\n')



    def write_init(self, output_file='BFMcoupler_VARDIAGinitializ.h'):

        # First part of the code
        header = '''C $Header:/MITgcm/pkg/BFMcoupler/BFMcoupler_VARDIAGinitializ.h,v 1.01
C initialization of local variables of diagnostic in BFMcoupler_calc_tendency.f
C 
      DO j=1-OLy,sNy+OLy
      DO i=1-OLx,sNx+OLx
          do k=1,Nr
'''

        footer = '''          enddo
      ENDDO
      ENDDO'''

        # write to file
        with open(output_file, 'w') as ofile:

            ofile.write(header)

            # loop on diagnostic variables
            for var in self._diag_vars:
                ofile.write('              dia' + var[0] + '(i,j,k)=0. _d 0\n')

            # loop on MITgcm-BFM-only diagnostic variables
            for var in self._diag_vars_MIT:
                ofile.write('              dia' + var[1][0] + '(i,j,k)=0. _d 0\n')

            ofile.write(footer)



    def write_copy_from_d(self, output_file='BFMcoupler_VARDIAGcopy_fromD.h'):

        # First part of the code
        header = '''C $Header:/MITgcm/pkg/BFMcoupler/BFMcoupler_VARDIAGcopy_fromD.h,v 1.01
C copy OUTPUT_ECOLOGY diagnostics from vector d to local variables of diagnostic in BFMcoupler_calc_tendency.f
C '''

        # write to file
        with open(output_file, 'w') as ofile:
            
            ofile.write(header)

            # loop on diagnostic variables
            for var in self._diag_vars:
                ofile.write('              dia' + var[0] + '(i,j,1:kBot(i,j))=d(' + str(self._diag_vars.index(var)+1) + ',1:kBot(i,j))\n')

            # loop on MITgcm-BFM-only diagnostic variables
            for var in self._diag_vars_MIT:
                ofile.write('              dia' + var[1][0] + '(i,j,1:kBot(i,j))=er(1:kBot(i,j),' + str(var[0]) + ')\n')



    def write_fill_diags(self, output_file='BFMcoupler_VARDIAG_fill_diags.h'):

        # First part of the code
        header = '''C $Header:/MITgcm/pkg/BFMcoupler/BFMcoupler_VARDIAG_fill_diags.h,v 1.01
C fill the diagnostic memory using DIAGNOSTICS_FILL
'''

        with open(output_file, 'w') as ofile:

            ofile.write(header)

            # loop on diagnostic variables
            for var in self._diag_vars:
                varname=var[0] + " "*8
                ofile.write('        CALL DIAGNOSTICS_FILL(dia' + var[0] + ',\'' + varname[:8] + '\'\n')
                ofile.write('     &  ,0,Nr,2,bi,bj,myThid)\n')

            # loop on MITgcm-BFM-only diagnostic variables
            for var in self._diag_vars_MIT:
                ofile.write('        CALL DIAGNOSTICS_FILL(dia' + var[1][0] + ',\'' + var[1][0] + '\'\n')
                ofile.write('     &  ,0,Nr,2,bi,bj,myThid)\n')



    def write_data_ptracers(self, output_file='data.ptracers'):

        header = ' &PTRACERS_PARM01\n' + ' PTRACERS_numInUse=' + str(len(self._vars)) + ',\n PTRACERS_Iter0= 0,'

        with open(output_file, 'w') as ofile:

            ofile.write(header)

            for var in self._vars:

                idx = str(self._vars.index(var))

                ofile.write('# tracer ' + idx + ' - ' + var[0] + '\n')
                ofile.write(' PTRACERS_names(' + idx + ')=\'' + var[0] + '\',\n')
                ofile.write(' PTRACERS_long_names(' + idx + ')=\'' + var[0] + '\',\n')
                ofile.write(' PTRACERS_units(' + idx + ')=\'' + var[1] + '\',\n')
                
                if 'ADVscheme' in self._vars_properties:
                    ofile.write(' PTRACERS_advScheme(' + idx + ')=' + self._vars_properties['ADVscheme'] + ',\n')
                if 'diffKh' in self._vars_properties:
                    ofile.write(' PTRACERS_diffKh(' + idx + ')=' + self._vars_properties['diffKh'] + ',\n')
                if 'diffKr' in self._vars_properties:
                    ofile.write(' PTRACERS_diffKr(' + idx + ')=' + self._vars_properties['diffKr'] + ',\n')
                if 'useGMRedi' in self._vars_properties:
                    ofile.write(' PTRACERS_useGMRedi(' + idx + ')=' + self._vars_properties['useGMRedi'] + ',\n')
                if 'useKPP' in self._vars_properties:
                    ofile.write(' PTRACERS_useKPP(' + idx + ')=' + self._vars_properties['useKPP'] + ',\n')
                if 'initialFile' in self._vars_properties:
                    ofile.write(' PTRACERS_initialFile(' + idx + ')=\'' + self._vars_properties['initialFile'] + var[0] + '.bin\'' + ',\n')
                if 'EvPrRn' in self._vars_properties:
                    ofile.write('#PTRACERS_EvPrRn(' + idx + ')=' + self._vars_properties['EvPrRn'] + ',\n')

                ofile.write('#PTRACERS_ref(1:Nr,' + idx + ')=Null,Null,Null ... ,\n')
                ofile.write(' &END\n')



    def write_data_diagnostic(self, output_file='data.diagnostic'):

        default_frequency = '432000'
        default_time_phase = '0'

        header = '''# Diagnostic Package Choices
#--------------------
#  dumpAtLast (logical): always write output at the end of simulation (default=F)
#  diag_mnc   (logical): write to NetCDF files (default=useMNC)
#--for each output-stream:
#  fileName(n) : prefix of the output file name (max 80c long) for outp.stream n
#  frequency(n):< 0 : write snap-shot output every |frequency| seconds
#               > 0 : write time-average output every frequency seconds
#  timePhase(n)     : write at time = timePhase + multiple of |frequency|
#    averagingFreq  : frequency (in s) for periodic averaging interval
#    averagingPhase : phase     (in s) for periodic averaging interval
#    repeatCycle    : number of averaging intervals in 1 cycle
#  levels(:,n) : list of levels to write to file (Notes: declared as REAL)
#                when this entry is missing, select all common levels of this list
#  fields(:,n) : list of selected diagnostics fields (8.c) in outp.stream n
#                (see "available_diagnostics.log" file for the full list of diags)
#  missing_value(n) : missing value for real-type fields in output file "n"
#  fileFlags(n)     : specific code (8c string) for output file "n"
#--------------------
 &DIAGNOSTICS_LIST
# diag_mnc     = .FALSE.,
# dumpAtLast   = .TRUE.,
'''

        footer = ''' &END
  
# Parameter for Diagnostics of per level statistics:
#-----------------
# for each output-stream:
#  stat_fname(n) : prefix of the output file name (only 8.c long) for outp.stream n
#  stat_freq(n):< 0 : write snap-shot output every |stat_freq| seconds
#               > 0 : write time-average output every stat_freq seconds
#  stat_phase(n)    : write at time = stat_phase + multiple of |stat_freq|
#  stat_region(:,n) : list of "regions" (default: 1 region only=global)
#  stat_fields(:,n) : list of diagnostics fields (8.c) (see "available_diagnostics.log"
#                 file for the list of all available diag. in th
#-----------------
 &DIAG_STATIS_PARMS
 &END
#--
'''

        with open(output_file, 'w') as ofile:

            ofile.write(header)

            # loop on diagnostic variables
            for var in self._diag_vars:
                ofile.write(' fields(1,' + str(self._diag_vars.index(var)) + ')  = \'' + var[0] + '\',\n')
                ofile.write(' fileName(' + str(self._diag_vars.index(var)) + ')  = \'' + var[0] + '\',\n')
                ofile.write(' frequency(' + str(self._diag_vars.index(var)) + ')  = ' + default_frequency + ',\n')
                ofile.write(' timePhase(' + str(self._diag_vars.index(var)) + ')  = ' + default_time_phase + ',\n')
                ofile.write('#\n')

            # loop on MITgcm-BFM-only diagnostic variables
            for var in self._diag_vars_MIT:
                ofile.write(' fields(1,' + str(var[0]) + ')  = \'' + var[1][0] + '\',\n')
                ofile.write(' fileName(' + str(var[0]) + ')  = \'' + var[1][0] + '\',\n')
                ofile.write(' frequency(' + str(var[0]) + ')  = ' + default_frequency + ',\n')
                ofile.write(' timePhase(' + str(var[0]) + ')  = ' + default_time_phase + ',\n')
                ofile.write('#\n')

            ofile.write(footer)



if __name__ == '__main__':

    # main object istantiation
    my_BFM_vars = BFM_vars(args.inputfile)

    # add missing MITgcm-BFM variables
    #my_BFM_vars.add_diag_var_MIT('wspeed', 'm/s', 9)
    #my_BFM_vars.add_diag_var_MIT('PCO2atm', 'kg/m3', 5)

    # define additional tracers properties
    my_BFM_vars.set_vars_property('ADVscheme', '33')
    my_BFM_vars.set_vars_property('diffKh', '400')
    my_BFM_vars.set_vars_property('diffKr', '1.E-4')
    my_BFM_vars.set_vars_property('useGMRedi', '.FALSE.')
    my_BFM_vars.set_vars_property('useKPP', '.FALSE.')
    my_BFM_vars.set_vars_property('initialFile', '../input/input_binaries/IC_SAD_')
    my_BFM_vars.set_vars_property('EvPrRn', '0.0')

    # debug
    #my_BFM_vars.flush()

    # output files
    if args.type=='code':
        my_BFM_vars.write_diagnostic_init(OUTDIR+'BFMcoupler_diagnostics_init.F')
        my_BFM_vars.write_local(OUTDIR+'BFMcoupler_VARDIAGlocal.h')
        my_BFM_vars.write_init(OUTDIR+'BFMcoupler_VARDIAGinitializ.h')
        my_BFM_vars.write_copy_from_d(OUTDIR+'BFMcoupler_VARDIAGcopy_fromD.h')
        my_BFM_vars.write_fill_diags(OUTDIR+'BFMcoupler_VARDIAG_fill_diags.h')
    if args.type=='namelist':
        my_BFM_vars.write_data_ptracers(OUTDIR + 'data.ptracers')
        my_BFM_vars.write_data_diagnostic(OUTDIR + 'data.diagnostic')
