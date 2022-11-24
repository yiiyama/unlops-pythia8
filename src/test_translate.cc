#include "Pythia8/MergingHooks.h"
#include <iostream>

int
main(int argc, char* argv[])
{
  Pythia8::HardProcess proc;
  proc.translateProcessString(argv[1]);

  std::cout << "nQ " << proc.nQuarksOut() << std::endl;
  std::cout << "nL " << proc.nLeptonOut() << std::endl;
  std::cout << "nB " << proc.nBosonsOut() << std::endl;

  return 0;
}
