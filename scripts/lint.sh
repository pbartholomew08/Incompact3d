#!/bin/bash
#
#        FILE: lint.sh
#      AUTHOR: Paul Bartholomew <ptb08@ic.ac.uk>
# DESCRIPTION: Runs a 'lint' build of Incompact3D, with no optimisation but maximum warnings.
#              Output in forms of warnings is written to lint.txt
#              Note this script should be run from the Incompact3D root directory.
#

# Clean, run lint build and cleanup
make clean
make lint 2> lint.txt
make clean
