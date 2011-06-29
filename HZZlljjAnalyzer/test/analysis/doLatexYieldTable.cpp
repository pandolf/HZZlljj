#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include "TFile.h"
#include "TH1D.h"


struct HistoFiles {

  TFile* file_signal250;
  TFile* file_signal300;
  TFile* file_signal350;
  TFile* file_signal400;
  TFile* file_signal450;
  TFile* file_signal500;
  TFile* file_ZJets;
  TFile* file_Diboson;
  TFile* file_TT;
  TFile* file_ZCC;
  TFile* file_ZBB;

};





void writeTableFile( HistoFiles files, const std::string& selectionType, const std::string& btagcat );
std::string getSingleLine( int mass, const std::string& btagcat, HistoFiles files );



int main( int argc, char* argv[] ) {

  if( argc != 2 ) {
    std::cout << "USAGE: ./doLatexYieldTable.cpp [selectionType]" << std::endl;
    return 1;
  }

  std::string selectionType(argv[1]); 
  

  HistoFiles files;

  std::string fileName;

  fileName = "HZZlljjRM_SMHiggsToZZTo2L2Q_M-250_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2_" + selectionType + "_ALL.root";
  TFile* file_signal250 = TFile::Open(fileName.c_str());
  fileName = "HZZlljjRM_SMHiggsToZZTo2L2Q_M-300_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2_" + selectionType + "_ALL.root";
  TFile* file_signal300 = TFile::Open(fileName.c_str());
  fileName = "HZZlljjRM_SMHiggsToZZTo2L2Q_M-350_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2_" + selectionType + "_ALL.root";
  TFile* file_signal350 = TFile::Open(fileName.c_str());
  fileName = "HZZlljjRM_SMHiggsToZZTo2L2Q_M-400_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2_" + selectionType + "_ALL.root";
  TFile* file_signal400 = TFile::Open(fileName.c_str());
  fileName = "HZZlljjRM_SMHiggsToZZTo2L2Q_M-450_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2_" + selectionType + "_ALL.root";
  TFile* file_signal450 = TFile::Open(fileName.c_str());
  fileName = "HZZlljjRM_SMHiggsToZZTo2L2Q_M-500_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2_" + selectionType + "_ALL.root";
  TFile* file_signal500 = TFile::Open(fileName.c_str());
  fileName = "HZZlljjRM_ZJets_alpgen_TuneZ2_Spring11_v2_" + selectionType + "_ALL.root";
  TFile* file_ZJets = TFile::Open(fileName.c_str());
  fileName = "HZZlljjRM_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11_v2_" + selectionType + "_ALL.root";
  TFile* file_Diboson = TFile::Open(fileName.c_str());
  fileName = "HZZlljjRM_TT_TW_TuneZ2_7TeV-pythia6-tauola_Spring11_v2_" + selectionType + "_ALL.root";
  TFile* file_TT = TFile::Open(fileName.c_str());
  fileName = "HZZlljjRM_ZCC_alpgen_TuneZ2_Spring11_v2_" + selectionType + "_ALL.root";
  TFile* file_ZCC = TFile::Open(fileName.c_str());
  fileName = "HZZlljjRM_ZBB_alpgen_TuneZ2_Spring11_v2_" + selectionType + "_ALL.root";
  TFile* file_ZBB = TFile::Open(fileName.c_str());

  files.file_signal250 = file_signal250;
  files.file_signal300 = file_signal300;
  files.file_signal350 = file_signal350;
  files.file_signal400 = file_signal400;
  files.file_signal450 = file_signal450;
  files.file_signal500 = file_signal500;
  files.file_ZJets = file_ZJets;
  files.file_Diboson = file_Diboson;
  files.file_TT = file_TT;
  files.file_ZCC = file_ZCC;
  files.file_ZBB = file_ZBB;

  writeTableFile( files, selectionType, "0btag" );
  writeTableFile( files, selectionType, "1btag" );
  writeTableFile( files, selectionType, "2btag" );

  return 0;

}



void writeTableFile( HistoFiles files, const std::string& selectionType, const std::string& btagcat ) {


  std::string fileName = "yields_latex_" + selectionType + "_" + btagcat + ".txt";

  std::ofstream ofs(fileName.c_str());

  ofs << "\\begin{center}" << std::endl;
  ofs << "\\begin{tabular}{cccccccc}" << std::endl;
  ofs << "\\hline" << std::endl;
  ofs << "Mass [GeV] & Signal & Z+Jets & Z+cc & Z+bb & Diboson & tt/tW & Total Background \\\\" << std::endl;
  ofs << "\\hline" << std::endl;

  ofs << getSingleLine( 250, btagcat, files ) << std::endl;
  ofs << getSingleLine( 300, btagcat, files ) << std::endl;
  ofs << getSingleLine( 350, btagcat, files ) << std::endl;
  ofs << getSingleLine( 400, btagcat, files ) << std::endl;
  ofs << getSingleLine( 450, btagcat, files ) << std::endl;
  ofs << getSingleLine( 500, btagcat, files ) << std::endl;

  ofs << "\\hline" << std::endl;
  ofs << "\\end{tabular}" << std::endl;

  ofs.close();

}


std::string getSingleLine( int mass, const std::string& btagcat, HistoFiles files ) {

  char yieldHistoName[300];
  sprintf( yieldHistoName, "nEvents_fb_%s_%d", btagcat.c_str(), mass );

  char effHistoName[300];
  sprintf( effHistoName, "eff_%s_%d", btagcat.c_str(), mass );

  TFile* signalFile;
  if( mass==250 ) signalFile = files.file_signal250;
  else if( mass==300 ) signalFile = files.file_signal300;
  else if( mass==350 ) signalFile = files.file_signal350;
  else if( mass==400 ) signalFile = files.file_signal400;
  else if( mass==450 ) signalFile = files.file_signal450;
  else if( mass==500 ) signalFile = files.file_signal500;
  else {
    std::cout << "Unkown mass: " << mass << ". Exiting." << std::endl;
    exit(11);
  }


  TH1D* h1_yield_signal = (TH1D*)signalFile->Get(yieldHistoName);
  TH1D* h1_eff_signal = (TH1D*)signalFile->Get(effHistoName);

  TH1D* h1_yield_ZJets = (TH1D*)files.file_ZJets->Get(yieldHistoName);
  TH1D* h1_yield_ZCC = (TH1D*)files.file_ZCC->Get(yieldHistoName);
  TH1D* h1_yield_ZBB = (TH1D*)files.file_ZBB->Get(yieldHistoName);
  TH1D* h1_yield_Diboson = (TH1D*)files.file_Diboson->Get(yieldHistoName);
  TH1D* h1_yield_TT = (TH1D*)files.file_TT->Get(yieldHistoName);


  float yield_signal = h1_yield_signal->GetBinContent(1);
  float eff_signal = h1_eff_signal->GetBinContent(1);

  float yield_ZJets   = h1_yield_ZJets->GetBinContent(1);
  float yield_ZCC     = h1_yield_ZCC->GetBinContent(1);
  float yield_ZBB     = h1_yield_ZBB->GetBinContent(1);
  float yield_Diboson = h1_yield_Diboson->GetBinContent(1);
  float yield_TT      = h1_yield_TT->GetBinContent(1);

  float totalBG = yield_ZJets+yield_ZCC+yield_ZBB+yield_Diboson+yield_TT;

  std::setprecision(2);

  char returnLine[800];
  sprintf( returnLine, "%d & %.2g (%.2f\%%) & %.2g & %.2g & %.2g & %.2g & %.2g & %.2g \\\\", mass, yield_signal, 100.*eff_signal, yield_ZJets, yield_ZCC, yield_ZBB, yield_Diboson, yield_TT, totalBG);
  //sprintf( returnLine, "%d & %f (%f%%) & %f & %f & %f & %f & %f & %f \\\\", mass, yield_signal, eff_signal, yield_ZJets, yield_ZCC, yield_ZBB, yield_Diboson, yield_TT, totalBG);

  std::string returnLine_str(returnLine);

  return returnLine_str;

}

