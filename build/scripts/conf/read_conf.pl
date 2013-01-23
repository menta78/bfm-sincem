#!/usr/bin/perl -w

#Author: Esteban Gutierrez esteban.gutierrez@cmcc.it

use strict;
use warnings;

use Data::Dumper;
use Getopt::Std;

use read_memLayout; #reomve trailing white spaces
use read_namelist; # get the namelist values in the hash
use print_f90; # write the variables in the output file
use classes;


my ($nml_def, $nml_val, $proto_dir, $out_dir, @cpp_defs);

#fix values
my $VERBOSE = 0;
my $HELP    = 0;
my @PROTOS_NAME = qw(ModuleMem AllocateMem set_var_info_bfm init_var_bfm INCLUDE);
my @PROTOS_EXT  = qw(F90       F90         F90              F90          h      );

#structures to allocate the parameters read in memmory layour
my %lst_group = ();
my %lst_param = ();
my %lst_sta   = ();
my %lst_const = ();

sub usage(){
    print "usage: $0 {-D[cpp_def] -r [mem_layout] -n [namelist] -f [prototype_dir] -t [output_dir]} [-v]\n\n";
    print "This script generate .F90 and .h files using templates based on configuration files\n\n";
    print "MUST specify at least one these OPTIONS:\n";
    print "\t-D[cpp_def]          defines\n";
    print "\t-r [mem_layout]      memory layout\n";
    print "\t-n [namelist]        namelist\n";
    print "\t-f [prototype_dir]   input dir for prototype files\n";
    print "\t-t [output_dir]      output dir for generated files\n";
    print "alternative OPTIONS are:\n";
    print "\t-v                   verbose mode\n";
}

use Getopt::Long qw(:config bundling noignorecase); # for getopts compat
GetOptions(
    'r=s'  => \$nml_def,
    'n=s'  => \$nml_val,
    'f=s'  => \$proto_dir,
    't=s'  => \$out_dir,  
    'D=s@' => \@cpp_defs,
    'v'    => \$VERBOSE,
    'h'    => \$HELP,
    ) or &usage() && exit;
if ( $HELP ){ &usage(); exit; }
if ( !$nml_def || !$nml_val || !$proto_dir || !$out_dir || !@cpp_defs ){ &usage(); exit; }

#read memory layout file
read_memLayout($nml_def, \%lst_group, \%lst_param, \%lst_sta, \%lst_const, join(' ',@cpp_defs) );
#read_namelist($nml_val, \%lst_group, \%lst_param);
foreach my $idx ( 0 .. $#PROTOS_NAME ){
    my $name = $PROTOS_NAME[$idx];
    my $ext  = $PROTOS_EXT[$idx];
    if( $VERBOSE ){ print "---------- ${name} ini \n"; }
    print_f90("$proto_dir/${name}.proto", "$out_dir/${name}.${ext}", \%lst_group, \%lst_param, \%lst_sta, \%lst_const, $VERBOSE);
    if( $VERBOSE ){ print "---------- ${name} end \n\n"; }
}

#my @const = sort keys %lst_const; print "Constituents: @const\n";
#foreach my $key (sort keys %lst_group) { $lst_group{$key}->print(); }
#foreach my $key (sort keys %lst_param) { $lst_param{$key}->print(); }
#foreach my $key (sort keys %lst_sta) { print "STADISTICS => $key have $lst_sta{$key} elements\n"; }
#foreach my $key (sort keys %lst_param) { if( $lst_param{$key}->getType() eq 'flux' ){ $lst_param{$key}->print(); } }
#foreach my $key (sort keys %lst_param) { if( $lst_param{$key}->getInclude() ){ $lst_param{$key}->print(); } }
#print Dumper(\%lst_sta) , "\n";
#print Dumper(\%lst_const) , "\n";
