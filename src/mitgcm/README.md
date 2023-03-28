# BFMCOUPLER
MITgcm-BFM coupler package

### Intended workflow:

- compile bfmv5;
- modify `$BFM_ROOT_DIR/include/BFM_var_list.h` in order to fit fortran77 standard;
- launch `bfm_config_gen.py` script, in order to generate missing headers and source files accordingly to `namelist.passivetrc`;
- move generated files into the proper MIT code directory (as usual);
- compile and run MIT.
