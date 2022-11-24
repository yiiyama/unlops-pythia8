CXX = g++

CFLAGS_STD = -std=c++1y
CFLAGS_ROOT = -I$(shell root-config --incdir)
CFLAGS_PYTHIA8 = -I/afs/cern.ch/work/y/yiiyama/ext/include
CFLAGS_HEPMC2 = -I/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/hepmc/2.06.07-oenich/include
CFLAGS_FASTJET = $(shell fastjet-config --cxxflags)
CFLAGS = $(CFLAGS_ROOT) $(CFLAGS_PYTHIA8) $(CFLAGS_HEPMC2) $(CFLAGS_FASTJET)

LFLAGS_ROOT = $(shell root-config --libs)
LFLAGS_PYTHIA8 = -L/afs/cern.ch/work/y/yiiyama/ext/lib -lpythia8
LFLAGS_HEPMC2 = -L/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/hepmc/2.06.07-oenich/lib -lHepMC -lHepMCfio
LFLAGS_FASTJET = $(shell fastjet-config --libs)
LFLAGS_TBB = -L/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/tbb/2017_20161004oss/lib -ltbb
LFLAGS = $(LFLAGS_ROOT) $(LFLAGS_PYTHIA8) $(LFLAGS_HEPMC2) $(LFLAGS_FASTJET) $(LFLAGS_TBB)

main: bin/run_unlops

bin/run_unlops: src/run_unlops.cc
	mkdir -p bin
	$(CXX) $(CFLAGS) -o $@ $^ $(LFLAGS)

bin/debug_unlops: src/debug_unlops.cc
	mkdir -p bin
	$(CXX) $(CFLAGS) -o $@ $^ $(LFLAGS)

bin/compute_xsec: src/compute_xsec.cc
	mkdir -p bin
	$(CXX) $(CFLAGS_STD) $(CFLAGS_ROOT) $(CFLAGS_PYTHIA8) -o $@ $^ $(LFLAGS_ROOT) $(LFLAGS_PYTHIA8)

bin/debug_wgt: src/debug_wgt.cc
	mkdir -p bin
	$(CXX) $(CFLAGS_STD) $(CFLAGS_ROOT) $(CFLAGS_PYTHIA8) -o $@ $^ $(LFLAGS_ROOT) $(LFLAGS_PYTHIA8)

bin/run_fxfx: src/run_fxfx.cc
	mkdir -p bin
	$(CXX) $(CFLAGS) -o $@ $^ $(LFLAGS)

bin/test_translate: src/test_translate.cc
	mkdir -p bin
	$(CXX) $(CFLAGS_STD) $(CFLAGS_PYTHIA8) -o $@ $^ $(LFLAGS_PYTHIA8)

clean:
	rm -f bin/*
