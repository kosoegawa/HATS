#!/bin/bash

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# last modified and documented on February 17 2026
# This script was designed to run all run scripts


./runHlaA.pl
./runHlaB.pl
./runHlaC.pl

./runDRB1.pl
./runDRB3.pl
./runDRB4.pl
./runDRB5.pl

./runDQB1.pl
./runDQA1.pl

./runDPB1.pl
./runDPA1.pl

./combine.pl

./runBw46.pl

exit 0
