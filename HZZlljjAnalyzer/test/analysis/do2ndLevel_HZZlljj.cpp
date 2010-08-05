#include "Ntp1Analyzer_HZZlljj.h"


int main() {

  Ntp1Analyzer_HZZlljj* na = new Ntp1Analyzer_HZZlljj("H130");

  na->LoadInput();
  na->Loop();

  delete na;

  return 0;

}
