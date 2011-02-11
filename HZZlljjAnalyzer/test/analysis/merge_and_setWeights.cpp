#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TRegexp.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <cstdlib>


TChain* tree = 0;
std::string analysisType_;
std::string flags_;

struct EventsAndLumi {
  int nTotalEvents;
  float totalLumi;
};



EventsAndLumi addInput( const std::string& dataset );
float getWeight( const std::string& dataset, int nEvents );


int main( int argc, char* argv[] ) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./merge_and_setWeights [dataset] [analysisType=\"HZZlljj\"] [flags=\"\"]" << std::endl;
    exit(917);
  }

  std::string dataset = argv[1];

  analysisType_ = "HZZlljj";
  if( argc>=3 ) {
    std::string analysisType_str(argv[2]);
    analysisType_ = analysisType_str;
  }

  flags_ = "";
  if( argc==4 ) {
    std::string flags_str(argv[3]);
    flags_ = flags_str;
  }


  tree = new TChain("reducedTree");

  EventsAndLumi evlu;
  evlu = addInput( dataset );

  std::cout << std::endl << "-> Finished adding. Total entries: " << tree->GetEntries() << std::endl;

  float weight = getWeight( dataset, evlu.nTotalEvents );

  // and now set the weights
  tree->SetBranchStatus( "eventWeight", 0 );
  
  std::string outfilename = analysisType_ + "_2ndLevelTreeW_"+dataset;
  if( flags_!="" ) outfilename += "_" + flags_;
  outfilename += +".root";
  TFile* outfile = new TFile(outfilename.c_str(), "recreate");
  outfile->cd();

  TH1F* h1_lumi = new TH1F("lumi", "", 1, 0., 1.);
  h1_lumi->SetBinContent(1, evlu.totalLumi);
  TH1F* h1_nCounter = new TH1F("nCounter", "", 1, 0., 1.);
  h1_nCounter->SetBinContent(1, evlu.nTotalEvents);

  TTree* newTree = tree->CloneTree(0);
  Float_t newWeight;
  newTree->Branch( "eventWeight", &newWeight, "newWeight/F" );

  int nentries = tree->GetEntries();
  for( unsigned ientry = 0; ientry<nentries; ++ientry ) {

    tree->GetEntry(ientry);

    if( (ientry % 100000) ==0 ) std::cout << "Entry: " << ientry << " /" << nentries << std::endl;

    newWeight = weight;

    if( dataset=="MU_Run2010B_PromptReco_v2_runs146240_146733" ) newWeight = 0.5;

    newTree->Fill();

  } //for entries

  h1_lumi->Write();
  h1_nCounter->Write();
  newTree->Write();
  outfile->Write();
  outfile->Close();

  return 0;

}


EventsAndLumi addInput( const std::string& dataset ) {

  std::string infileName = "files_"+analysisType_+"_2ndLevel_" + dataset;
  if( flags_!="" ) infileName += "_" + flags_;
  infileName += ".txt";
  TH1F* h1_lumi;
  TH1F* h1_nCounter;

  int totalEvents = 0;
  float totalLumi = 0.;

  //open from file.txt:
  FILE* iff = fopen(infileName.c_str(),"r");
  if(iff == 0) {
    std::cout << "cannot open input file '" << infileName << "' ... adding single file." << std::endl;
    infileName = analysisType_+"_2ndLevelTree_" + dataset;
    if( flags_!="" ) infileName += "_" + flags_;
    infileName += ".root";
    std::string treeName = infileName +"/reducedTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    TFile* infile = TFile::Open(infileName.c_str(), "READ");
    h1_nCounter = (TH1F*)infile->Get("nCounter");
    if( h1_nCounter!=0 ) {
      totalEvents += h1_nCounter->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no nCounter information. Skipping." << std::endl;
    }
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
      std::cout << "\tTotal lumi: " << totalLumi << " ub-1";
    } else {
      //std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping.";
    }
    std::cout << std::endl;
    infile->Close();

  } else { //if file is good:

    char singleLine[500];
    std::cout << "-> Correctly opened file: '" << infileName << "'." << std::endl;

    while( fscanf(iff, "%s", singleLine) !=EOF ) {

      std::string rootfilename(singleLine);
      std::string treename = rootfilename + "/reducedTree";
      std::cout << "-> Added " << treename;
      tree->Add(treename.c_str());
      TFile* infile = TFile::Open(rootfilename.c_str(), "READ");
      h1_nCounter = (TH1F*)infile->Get("nCounter");
      if( h1_nCounter!=0 ) {
        totalEvents += h1_nCounter->GetBinContent(1);
      } else {
        std::cout << std::endl << " WARNING! File '" << rootfilename << "' has no nCounter information. Skipping." << std::endl;
      }
      h1_lumi = (TH1F*)infile->Get("lumi");
      if( h1_lumi!=0 ) {
        totalLumi += h1_lumi->GetBinContent(1);
        std::cout << "\tTotal lumi: " << totalLumi << " ub-1";
      } else {
        //std::cout << std::endl << " WARNING! File '" << rootfilename << "' has no lumi information. Skipping." << std::endl;
      }
      std::cout << std::endl;
      infile->Close();

    }
    fclose(iff);

  }

  EventsAndLumi evlu;
  evlu.nTotalEvents = totalEvents;
  evlu.totalLumi = totalLumi;

  return evlu;

} //addinput


float getWeight( const std::string& dataset, int nEvents ) {

  float xSection = -1.;

  // all cross sections in pb-1:
  if( dataset=="ZJets_madgraph" ) {
    xSection = 3048.; //NNLO see https://twiki.cern.ch/twiki/pub/CMS/GeneratorMain/ShortXsec.pdf
  } else if( dataset=="Z0Jets_Pt0to100-alpgen_Spring10" ) {
    xSection = 2350.*0.853 ; // sigma x filter efficiency taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/ProductionReProcessingSpring10#ALPGEN
  } else if( dataset=="Z1Jets_Pt0to100-alpgen_Spring10" ) {
    xSection = 870.*0.447 ; // sigma x filter efficiency
  } else if( dataset=="Z1Jets_Pt100to300-alpgen_Spring10" ) {
    xSection = 19.3*0.504; // sigma x filter efficiency
  } else if( dataset=="Z1Jets_Pt300to800-alpgen_Spring10" ) {
    xSection = 0.226*0.352; // sigma x filter efficiency
  } else if( dataset=="Z1Jets_Pt800to1600-alpgen_Spring10" ) {
    xSection = 0.000528*0.266; // sigma x filter efficiency
  } else if( dataset=="Z2Jets_Pt0to100-alpgen_Spring10" ) {
    xSection = 372.*0.26; // sigma x filter efficiency
  } else if( dataset=="Z2Jets_Pt100to300-alpgen_Spring10" ) {
    xSection = 27.0*0.315; // sigma x filter efficiency
  } else if( dataset=="Z2Jets_Pt300to800-alpgen_Spring10" ) {
    xSection = 0.457*0.244; // sigma x filter efficiency
  } else if( dataset=="Z2Jets_Pt800to1600-alpgen_Spring10" ) {
    xSection = 0.00132*0.203; // sigma x filter efficiency
  } else if( dataset=="Z3Jets_Pt0to100-alpgen_Spring10" ) {
    xSection = 140.*0.157; // sigma x filter efficiency
  } else if( dataset=="Z3Jets_Pt100to300-alpgen_Spring10" ) {
    xSection = 20.3*0.189; // sigma x filter efficiency
  } else if( dataset=="Z3Jets_Pt300to800-alpgen_Spring10" ) {
    xSection = 0.465*0.162; // sigma x filter efficiency
  } else if( dataset=="Z3Jets_Pt800to1600-alpgen_Spring10" ) {
    xSection = 0.00152*0.149; // sigma x filter efficiency
  } else if( dataset=="Z4Jets_Pt0to100-alpgen_Spring10" ) {
    xSection = 46.1*0.0939; // sigma x filter efficiency
  } else if( dataset=="Z4Jets_Pt100to300-alpgen_Spring10" ) {
    xSection = 10.7*0.115; // sigma x filter efficiency
  } else if( dataset=="Z4Jets_Pt300to800-alpgen_Spring10" ) {
    xSection = 0.319*0.104; // sigma x filter efficiency
  } else if( dataset=="Z4Jets_Pt800to1600-alpgen_Spring10" ) {
    xSection = 0.0011*0.106; // sigma x filter efficiency
  } else if( dataset=="Z5Jets_Pt0to100-alpgen_Spring10" ) {
    xSection = 13.9*0.0727; // sigma x filter efficiency
  } else if( dataset=="Z5Jets_Pt100to300-alpgen_Spring10" ) {
    xSection = 4.42*0.0956; // sigma x filter efficiency
  } else if( dataset=="Z5Jets_Pt300to800-alpgen_Spring10" ) {
    xSection = 0.164*0.103; // sigma x filter efficiency
  } else if( dataset=="Z5Jets_Pt800to1600-alpgen_Spring10" ) {
    xSection = 0.000588*0.109; // sigma x filter efficiency
  } else if( dataset=="Z0Jets_TuneD6T" || dataset=="Z0Jets_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 1929. ; // cross sections taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/ProductionFall2010#ALPGEN
  } else if( dataset=="Z1Jets_ptZ-0to100_TuneD6T" || dataset=="Z1Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 380.8; 
  } else if( dataset=="Z1Jets_ptZ-100to300_TuneD6T" || dataset=="Z1Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 8.721;
  } else if( dataset=="Z1Jets_ptZ-300to800_TuneD6T" || dataset=="Z1Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 0.07386;
  } else if( dataset=="Z1Jets_ptZ-800to1600_TuneD6T" || dataset=="Z1Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 0.0001374;
  } else if( dataset=="Z2Jets_ptZ-0to100_TuneD6T" || dataset=="Z2Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 103.5;
  } else if( dataset=="Z2Jets_ptZ-100to300_TuneD6T" || dataset=="Z2Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 8.534;
  } else if( dataset=="Z2Jets_ptZ-300to800_TuneD6T" || dataset=="Z2Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 0.1151;
  } else if( dataset=="Z2Jets_ptZ-800to1600_TuneD6T" || dataset=="Z2Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 0.0003023;
  } else if( dataset=="Z3Jets_ptZ-0to100_TuneD6T" || dataset=="Z3Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 22.89; 
  } else if( dataset=="Z3Jets_ptZ-100to300_TuneD6T" || dataset=="Z3Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 3.951;
  } else if( dataset=="Z3Jets_ptZ-300to800_TuneD6T" || dataset=="Z3Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 0.08344;
  } else if( dataset=="Z3Jets_ptZ-800to1600_TuneD6T" || dataset=="Z3Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 0.0002480;
  } else if( dataset=="Z4Jets_ptZ-0to100_TuneD6T" || dataset=="Z4Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 4.619;
  } else if( dataset=="Z4Jets_ptZ-100to300_TuneD6T" || dataset=="Z4Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 1.298;
  } else if( dataset=="Z4Jets_ptZ-300to800_TuneD6T" || dataset=="Z4Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 0.03935;
  } else if( dataset=="Z4Jets_ptZ-800to1600_TuneD6T" || dataset=="Z4Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 0.0001394;
  } else if( dataset=="Z5Jets_ptZ-0to100_TuneD6T" || dataset=="Z5Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 1.135;
  } else if( dataset=="Z5Jets_ptZ-100to300_TuneD6T" || dataset=="Z5Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 0.4758;
  } else if( dataset=="Z5Jets_ptZ-300to800_TuneD6T" || dataset=="Z5Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 0.01946;
  } else if( dataset=="Z5Jets_ptZ-800to1600_TuneD6T" || dataset=="Z5Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola" ) {
    xSection = 0.00007195;
  } else if( dataset=="HZZ_qqll_gluonfusion_M130" ) {
    xSection = 25.560*0.03913*0.10097*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="HZZ_qqll_gluonfusion_M150" ) {
    xSection = 19.568*0.08234*0.10097*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="HZZ_qqll_gluonfusion_M200" ) {
    //xSection = 10.361*0.2537*0.10097*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
    xSection = 10.361*0.2537*0.10097*0.7*2.*40.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2 ENHANCE SIGMA BY 40!!!!!
  } else if( dataset=="HZZ_qqll_gluonfusion_M300" ) {
    xSection = 5.2728*0.3053*0.10097*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="HZZ_qqll_gluonfusion_M400" ) {
    xSection = 4.8236*0.2664*0.10097*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="HZZ_qqll_gluonfusion_M500" ) {
    xSection = 2.1914*0.2602*0.10097*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM200_2l2j_FASTSIM" || dataset=="SMHiggsToZZTo2L2Q_M-200_7TeV-jhu-pythia6" ) {
    xSection = (5.4218+0.64243+0.015972+0.056746+0.10830)*0.2537*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM250_2l2j_FASTSIM"  || dataset=="SMHiggsToZZTo2L2Q_M-250_7TeV-jhu-pythia6") {
    xSection = (3.4120+0.43490+0.0073067+0.023553+0.044954)*0.2951*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM300_2l2j_FASTSIM" || dataset=="SMHiggsToZZTo2L2Q_M-300_7TeV-jhu-pythia6" ) {
    xSection = (2.5001+0.30317+0.0038098+0.010992+0.021079)*0.3053*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM350_2l2j_FASTSIM" || dataset=="SMHiggsToZZTo2L2Q_M-350_7TeV-jhu-pythia6" ) {
    xSection = (2.4983+0.21635+0.0021911+0.0055876+0.010794)*0.3023*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM400_2l2j_FASTSIM" || dataset=="SMHiggsToZZTo2L2Q_M-400_7TeV-jhu-pythia6" ) {
    xSection = (2.0608+0.15739+0.0054837+0.00030308+0.0059045)*0.2664*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM450_2l2j_FASTSIM" || dataset=="SMHiggsToZZTo2L2Q_M-450_7TeV-jhu-pythia6" ) {
    xSection = (1.3445+0.11630+0.00087374+0.0017291+0.00033991)*0.2582*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM500_2l2j"|| dataset=="JHUgen_HiggsSM500_2l2j_FASTSIM" || dataset=="SMHiggsToZZTo2L2Q_M-500_7TeV-jhu-pythia6" ) {
    xSection = (0.83396+0.087068+0.00058401+0.0010272+0.0020377)*0.2602*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM550_2l2j_FASTSIM" || dataset=="SMHiggsToZZTo2L2Q_M-550_7TeV-jhu-pythia6" ) {
    xSection = (0.51418+0.065912+0.00039969+0.00063066+0.0012623)*0.2657*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM600_2l2j_FASTSIM" || dataset=="SMHiggsToZZTo2L2Q_M-600_7TeV-jhu-pythia6" ) {
    xSection = (0.31942+0.050376+0.00027833+0.00039795+0.00080323)*0.2724*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="TTbar_2l_Spring10" ) {
    xSection = 157.4*0.1080*0.1080*3.*3.; //NLO x BR(W->lnu) see https://twiki.cern.ch/twiki/pub/CMS/GeneratorMain/ShortXsec.pdf
  } else if( dataset=="TTJets_TuneD6T" || dataset=="TTJets_TuneZ2_7TeV-madgraph-tauola_Fall10" ) {
    xSection = 157.4; //NLO see https://twiki.cern.ch/twiki/pub/CMS/GeneratorMain/ShortXsec.pdf
  } else if( dataset=="ZZ_Spring10" || dataset=="ZZtoAnything_TuneZ2" || dataset=="ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10" ) {
    xSection = 5.9*1.3; //MCFM NLO see http://ceballos.web.cern.ch/ceballos/hwwlnln/cross_sections_backgrounds.txt plus factor 1.3 to account for glu-glu
  } else if( dataset=="WZtoAnything_TuneZ2" || dataset=="WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10" ) {
    xSection = 18.3; //MCFM NLO see http://ceballos.web.cern.ch/ceballos/hwwlnln/cross_sections_backgrounds.txt
  } else if( dataset=="WWtoAnything_TuneZ2" || dataset=="WWtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10" ) {
    xSection = 42.9; //MCFM NLO see http://ceballos.web.cern.ch/ceballos/hwwlnln/cross_sections_backgrounds.txt
  } else if( dataset=="Zmumu_Pythia" ) {
    xSection = 3048./3.; //NNLO see https://twiki.cern.ch/twiki/pub/CMS/GeneratorMain/ShortXsec.pdf
  } else if( dataset=="PhotonJet_Summer1036X_Pt5to15_pfakt5" ) {
    xSection = 4030000.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt15to20_pfakt5" ) {
    xSection = 114700.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt20to30_pfakt5" ) {
    xSection = 57180.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt30to50_pfakt5" ) { 
    xSection = 16520.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt50to80_pfakt5" ) {
    xSection = 2723.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt80to120_pfakt5" ) {
    xSection = 446.2;
  } else if( dataset=="PhotonJet_Summer1036X_Pt120to170_pfakt5" ) {
    xSection = 84.43;
  } else if( dataset=="PhotonJet_Summer1036X_Pt170to300_pfakt5" ) {
    xSection = 22.55;
  } else if( dataset=="PhotonJet_Summer1036X_Pt300to500_pfakt5" ) {
    xSection = 1.545;
  } else if( dataset=="PhotonJet_Summer1036X_Pt500toInf_pfakt5" ) {
    xSection = 0.0923;
  } else if( dataset=="QCD_Pt_120to170_TuneZ2_7TeV_pythia6" ) {
    xSection = 1.151e+05;
  } else if( dataset=="QCD_Pt_170to300_TuneZ2_7TeV_pythia6" ) {
    xSection = 2.426e+04;
  } else {
    std::cout << std::endl;
    std::cout << "-> WARNING!! Dataset: '" << dataset << "' not present in database. Cross section unknown." << std::endl;
    std::cout << "-> Will set unitary weights." << std::endl;
    return 1.;
  }

  TString dataset_tstr(dataset.c_str());
  TRegexp re("alpgen");
  TRegexp re2("Z");
  TRegexp re3("Jets");
  TRegexp re4("TuneD6T");
  if( dataset_tstr.Contains(re)||(dataset_tstr.Contains(re2)&&dataset_tstr.Contains(re3)&&dataset_tstr.Contains(re4)) ) {
    std::cout << "-> Scaling LO alpgen cross-section to NNLO." << std::endl;
    xSection*=1.31; // K factor
  }


  std::cout << "-> Total Events Analyzed: " << nEvents << ". " << std::endl;
  std::cout << "-> Dataset cross-section: " << xSection << " pb" << std::endl;
  float weight = xSection/((float)nEvents);
  std::cout << "=> Event weight: " << weight << std::endl;

  return weight;

}
