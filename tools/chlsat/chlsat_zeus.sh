#!/bin/bash
bashsrcdir=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
ulimit -s unlimited
. $bashsrcdir/set_env_nemo4.sh
$bashsrcdir/chlsat.x
