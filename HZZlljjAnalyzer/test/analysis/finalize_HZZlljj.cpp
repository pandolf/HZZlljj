
#include "Ntp1Finalizer_HZZlljj.h"
#include "TMath.h"
#include <iostream>



double delta_phi(double phi1, double phi2) {

  double dphi = fabs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
}


float delta_phi(float phi1, float phi2) {

  float dphi = fabs(phi1 - phi2);
  float sgn = (phi1 >= phi2 ? +1. : -1.);
  return sgn * (dphi <= TMath::Pi() ? dphi : TMath::TwoPi() - dphi);
}





int main( int argc, char* argv[] ) {

  if( argc!=2 && argc!=3 ) {
    std::cout << "USAGE: ./finalize_HZZlljj [dataset] [leptType=\"ALL\"]" <<std::endl;
    return 13;
  }


  std::string dataset(argv[1]);

  std::string leptType="ALL";
  if( argc==3 ) {
    std::string leptType_str(argv[2]);
    leptType = leptType_str;
  }



  Ntp1Finalizer_HZZlljj* nf = new Ntp1Finalizer_HZZlljj( dataset, leptType );



  if( dataset=="DATA_EG_37X" ) {

    nf->addFile( "EG_Run2010A_Jul15thReReco_v1" );
    nf->addFile( "EG_Run2010A_Jul26thReReco_v1" );

  } else if( dataset=="Electron_Run2010B" ) {

    nf->addFile( "Electron_Run2010B_PromptReco_v2_runs146240_146733_prod2" );
    nf->addFile( "Electron_Run2010B_PromptReco_v2_runs146733_147111" );

  } else if( dataset=="Run2010B_runs146240_146733" ) {

    nf->addFile( "Electron_Run2010B_PromptReco_v2_runs146240_146733_prod2" );
    nf->addFile( "MU_Run2010B_PromptReco_v2_runs146240_146733" );

  } else if( dataset=="Mu_upto147589" ) {

    nf->addFile("Mu_Run2010B_PromptReco_v2_runs146734_147589");
    nf->addFile("Mu_Run2010A-Sep17ReReco_v2_runs135821_144114");
    nf->addFile("Mu_Run2010B_PromptReco_v2_runs146240_146733");

  } else if( dataset=="ELEMUCombined" ) {

    nf->addFile("Mu_Run2010B_PromptReco_v2_runs146734_147589");
    nf->addFile("Mu_Run2010A-Sep17ReReco_v2_runs135821_144114");
    nf->addFile("Mu_Run2010B_PromptReco_v2_runs146240_146733");
    nf->addFile("EG_upto146724");

  } else if( dataset=="ZJets_alpgen" ) {

    nf->addFile( "Z0Jets_Pt0to100-alpgen_Spring10" );
    nf->addFile( "Z1Jets_Pt0to100-alpgen_Spring10" );
    nf->addFile( "Z1Jets_Pt100to300-alpgen_Spring10" );
    nf->addFile( "Z1Jets_Pt300to800-alpgen_Spring10" );
    nf->addFile( "Z1Jets_Pt800to1600-alpgen_Spring10" );
    nf->addFile( "Z2Jets_Pt0to100-alpgen_Spring10" );
    nf->addFile( "Z2Jets_Pt100to300-alpgen_Spring10" );
    nf->addFile( "Z2Jets_Pt300to800-alpgen_Spring10" );
    nf->addFile( "Z2Jets_Pt800to1600-alpgen_Spring10" );
    nf->addFile( "Z3Jets_Pt0to100-alpgen_Spring10" );
    nf->addFile( "Z3Jets_Pt100to300-alpgen_Spring10" );
    nf->addFile( "Z3Jets_Pt300to800-alpgen_Spring10" );
    nf->addFile( "Z3Jets_Pt800to1600-alpgen_Spring10" );
    nf->addFile( "Z4Jets_Pt0to100-alpgen_Spring10" );
    nf->addFile( "Z4Jets_Pt100to300-alpgen_Spring10" );
    nf->addFile( "Z4Jets_Pt300to800-alpgen_Spring10" );
    nf->addFile( "Z4Jets_Pt800to1600-alpgen_Spring10" );
    nf->addFile( "Z5Jets_Pt0to100-alpgen_Spring10" );
    nf->addFile( "Z5Jets_Pt100to300-alpgen_Spring10" );
    nf->addFile( "Z5Jets_Pt300to800-alpgen_Spring10" );
    nf->addFile( "Z5Jets_Pt800to1600-alpgen_Spring10" );

  } else if( dataset=="Z1Jets_alpgen_Spring10" ) {

    nf->addFile( "Z1Jets_Pt0to100-alpgen_Spring10" );
    nf->addFile( "Z1Jets_Pt100to300-alpgen_Spring10" );
    nf->addFile( "Z1Jets_Pt300to800-alpgen_Spring10" );
    nf->addFile( "Z1Jets_Pt800to1600-alpgen_Spring10" );

  } else if( dataset=="Z2Jets_alpgen_Spring10" ) {

    nf->addFile( "Z2Jets_Pt0to100-alpgen_Spring10" );
    nf->addFile( "Z2Jets_Pt100to300-alpgen_Spring10" );
    nf->addFile( "Z2Jets_Pt300to800-alpgen_Spring10" );
    nf->addFile( "Z2Jets_Pt800to1600-alpgen_Spring10" );

  } else if( dataset=="Z3Jets_alpgen_Spring10" ) {

    nf->addFile( "Z3Jets_Pt0to100-alpgen_Spring10" );
    nf->addFile( "Z3Jets_Pt100to300-alpgen_Spring10" );
    nf->addFile( "Z3Jets_Pt300to800-alpgen_Spring10" );
    nf->addFile( "Z3Jets_Pt800to1600-alpgen_Spring10" );

  } else if( dataset=="Z4Jets_alpgen_Spring10" ) {

    nf->addFile( "Z4Jets_Pt0to100-alpgen_Spring10" );
    nf->addFile( "Z4Jets_Pt100to300-alpgen_Spring10" );
    nf->addFile( "Z4Jets_Pt300to800-alpgen_Spring10" );
    nf->addFile( "Z4Jets_Pt800to1600-alpgen_Spring10" );

  } else if( dataset=="Z5Jets_alpgen_Spring10" ) {

    nf->addFile( "Z5Jets_Pt0to100-alpgen_Spring10" );
    nf->addFile( "Z5Jets_Pt100to300-alpgen_Spring10" );
    nf->addFile( "Z5Jets_Pt300to800-alpgen_Spring10" );
    nf->addFile( "Z5Jets_Pt800to1600-alpgen_Spring10" );

  } else if( dataset=="all" ) {


    if( leptType=="ELE" ) {
      argv[1] = new char[strlen("EG_upto146724")];
      strcpy( argv[1], "EG_upto146724");
      main( argc, argv );
    } else if( leptType=="MU") {
      argv[1] = new char[strlen("Mu_upto147589")];
      strcpy( argv[1], "Mu_upto147589");
      main( argc, argv );
      //finalize( "Mu_upto147589", leptType);
    } else if( leptType=="ALL" ) {
      argv[1] = new char[strlen("ELEMUCombined")];
      strcpy( argv[1], "ELEMUCombined");
      main( argc, argv );
    }


    argv[1] = new char[strlen("HZZ_qqll_gluonfusion_M130")];
    strcpy( argv[1], "HZZ_qqll_gluonfusion_M130");
    main( argc, argv );
    argv[1] = new char[strlen("HZZ_qqll_gluonfusion_M150")];
    strcpy( argv[1], "HZZ_qqll_gluonfusion_M150");
    main( argc, argv );
    argv[1] = new char[strlen("HZZ_qqll_gluonfusion_M200")];
    strcpy( argv[1], "HZZ_qqll_gluonfusion_M200");
    main( argc, argv );
    argv[1] = new char[strlen("HZZ_qqll_gluonfusion_M300")];
    strcpy( argv[1], "HZZ_qqll_gluonfusion_M300");
    main( argc, argv );
    argv[1] = new char[strlen("HZZ_qqll_gluonfusion_M400")];
    strcpy( argv[1], "HZZ_qqll_gluonfusion_M400");
    main( argc, argv );
    argv[1] = new char[strlen("HZZ_qqll_gluonfusion_M500")];
    strcpy( argv[1], "HZZ_qqll_gluonfusion_M500");
    main( argc, argv );
    argv[1] = new char[strlen("TTbar_2l_Spring10")];
    strcpy( argv[1], "TTbar_2l_Spring10");
    main( argc, argv );
    argv[1] = new char[strlen("ZZ_Spring10")];
    strcpy( argv[1], "ZZ_Spring10");
    main( argc, argv );
    argv[1] = new char[strlen("ZJets_madgraph")];
    strcpy( argv[1], "ZJets_madgraph");
    main( argc, argv );
    argv[1] = new char[strlen("ZJets_alpgen")];
    strcpy( argv[1], "ZJets_alpgen");
    main( argc, argv );
    //return 0;
    exit(1111);

  } else {
  
    nf->addFile( dataset );

  }

  std::cout << "-> Total integrated luminosity: " << nf->get_totalLumi() << " ub-1." << std::endl;

  nf->finalize();


  return 0;

}


