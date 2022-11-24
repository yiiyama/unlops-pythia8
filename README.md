# unlops-pythia8
Pythia8 programs to test and compare UNLOPS and FxFx generation

## Input LHE files

You will need to prepare the input LHE files for NLO 2j, 1j, 0j, LO 3j, 2j, and 1j processes by yourself.
To use the batch script `execute.sh`, all LHE files must contain the same number of events (assumed 40000
 in execute.sh). Lists of the LHE files should be compiled per process and be written into files with names
 Z_[N]LO_Inclusive_(0|1|2)j.list in some directory.

## Program compilation

Get the pythia source and have it built. Then edit the Makefile to suit your environment, then do `make`.

## Single run

Edit the `[CHANNELS]` block of cfg/run_unlops.cfg and do `bin/run_unlops cfg/run_unlops.cfg`.
