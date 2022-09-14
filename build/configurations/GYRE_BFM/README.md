# GYRE_BFM: the standard BFM-NEMO coupling

This version requires NEMO 4.2  and above.

It contains the following NEMO files in addition to the BFM files
- `namelist_cfg`
- `namelist_top_cfg`
- `iodef.xml`

Note that in case the NEMO code is upgraded, the previous files
must be manually changed accordingly

The run directory is created in `$BFMDIR_RUN/gyre_bfm` just like all other
BFM presets and the nemo executable is linked to there as well as
the NEMO reference configurations.

