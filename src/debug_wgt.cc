#include "Pythia8/Pythia.h"

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

class LHAupLHEFCustom : public Pythia8::LHAupLHEF {
public:
  LHAupLHEFCustom(Pythia8::Info* info, const char* filename) :
    LHAupLHEF(info, filename)
  {}

  bool useExternal() override { return true; }
  bool setEvent(int) override;

  double sumAbsW() const { return sumAbsW_; }

private:
  double sumAbsW_{0.};
};

bool
LHAupLHEFCustom::setEvent(int)
{
  bool success(setNewEventLHEF());

  sumAbsW_ += std::abs(xwgtupSave);

  if (!success)
    return false;
  else
    return setOldEventLHEF();
}

int
main(int argc, char* argv[])
{
  std::string configFile(argv[1]);
  std::string lheFile(argv[2]);
  unsigned order(std::atoi(argv[3]));
  unsigned nPartons(std::atoi(argv[4]));

  Pythia8::Pythia pythia("", false);

  pythia.readFile(configFile);

  pythia.settings.flag("PartonLevel:FSR", false);
  pythia.settings.flag("PartonLevel:ISR", false);
  pythia.settings.flag("HadronLevel:all", false);
  pythia.settings.flag("PartonLevel:MPI", false);

  // Switch on cross section estimation procedure.
  pythia.settings.flag("Merging:doXSectionEstimate", true);
  switch (order) {
  case 0:
    pythia.settings.flag("Merging:doUNLOPSTree", true);
    pythia.settings.flag("Merging:doUNLOPSLoop", false);
    break;
  case 1:
    pythia.settings.flag("Merging:doUNLOPSTree", false);
    pythia.settings.flag("Merging:doUNLOPSLoop", true);
    break;
  default:
    break;
  }

  pythia.settings.mode("Merging:nRequested", nPartons);

  pythia.settings.mode("Beams:frameType", 4);

  auto* lheInput = new LHAupLHEFCustom(&pythia.info, lheFile.c_str());
  pythia.setLHAupPtr(lheInput);

  pythia.init();

  double sumW(0.);
  double sumAbsW(0.);

  long long iEvent(0);
  long long nEvents(pythia.settings.mode("Main:numberOfEvents"));

  while (iEvent++ != nEvents) {
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile())
        break;
      continue;
    }

    sumW += pythia.info.weight();
    sumAbsW += std::abs(pythia.info.weight());
  }

  return 0;
}
