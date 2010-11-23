
#include "Ntp1Finalizer_ZGamma.h"
#include "TMath.h"
#include <iostream>






int main( int argc, char* argv[] ) {

  if( argc!=2 ) {
    std::cout << "USAGE: ./finalize_ZGamma [dataset]" <<std::endl;
    return 13;
  }


  std::string dataset(argv[1]);



  Ntp1Finalizer_ZGamma* nf = new Ntp1Finalizer_ZGamma( dataset );

  if( dataset=="PhotonJet_Summer1036X" ) {
    nf->addFile("PhotonJet_Summer1036X_Pt5to15_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt15to20_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt20to30_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt30to50_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt50to80_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt80to120_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt120to170_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt170to300_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt300to500_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt500toInf_pfakt5");
  } else {
    nf->addFile( dataset );
  }

 

  std::cout << "-> Total integrated luminosity: " << nf->get_totalLumi() << " ub-1." << std::endl;

  nf->finalize();

  delete nf;
  nf=0;

  return 0;

}


