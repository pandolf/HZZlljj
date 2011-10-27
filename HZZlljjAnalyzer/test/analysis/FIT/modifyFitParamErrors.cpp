#include <iostream>
#include <fstream>
#include <cstdlib>

#include <vector>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"

#include "SidebandFitter.h"




void modifyFitResults( const std::string& dataset, int nbtags, TFile* file, std::vector<std::string> varNames );


int main( int argc, char* argv[] ) {


  std::string dataset="LP11";
  if( argc>1 ) {
    std::string dataset_tmp(argv[1]);
    dataset = dataset_tmp;
  }


  std::string fitToysFileName = "fitParamErrors_" + dataset + ".root";

  std::string hadd_command = "hadd -f " + fitToysFileName + " /cmsrm/pc22_2/pandolf/FitToys_LP11/fitParamErrors_LP11_*";
  system(hadd_command.c_str());

  TFile* fitToysFile = TFile::Open(fitToysFileName.c_str());

  std::vector<std::string> varNames;
  varNames.push_back("alpha_rot");
  varNames.push_back("wdth_rot");


  modifyFitResults( dataset, 0, fitToysFile, varNames );
  modifyFitResults( dataset, 1, fitToysFile, varNames );
  modifyFitResults( dataset, 2, fitToysFile, varNames );

  return 0;

}




void modifyFitResults( const std::string& dataset, int nbtags, TFile* file, std::vector<std::string> varNames ) {


  // first: get errors
  std::vector<float> varErrors;
  for( unsigned ivar=0; ivar<varNames.size(); ++ivar ) {

    char histoName[500];
    sprintf( histoName, "%s_toys_%dbtag", varNames[ivar].c_str(), nbtags );

    TH1D* histo = (TH1D*)file->Get(histoName);

    TF1* gaussian = new TF1("gaussian", "gaus");
    histo->Fit(gaussian, "Q+");

    varErrors.push_back(gaussian->GetParameter(2));

  }

  SidebandFitter* sf = new SidebandFitter(dataset);

  std::string fitResultsFile_old = sf->get_fitResultsName( nbtags );
  std::string fitResultsFile_new = sf->get_fitResultsName( nbtags, "DATA_NEW" );

  std::ifstream ifs(fitResultsFile_old.c_str());
  std::ofstream ofs(fitResultsFile_new.c_str());

  ifs.clear();
  ifs.seekg(0);

  while( ifs.good() && !ifs.eof() ) {
  
    std::string varName;
    float value, error;

    ifs >> varName >> value >> error;

    bool foundVar=false;

    for( unsigned ivar=0; ivar<varNames.size() && !foundVar; ++ivar ) {
  
      if( varName==varNames[ivar] ) {
        ofs << varName << " " << value << " " << varErrors[ivar] << std::endl;
        foundVar = true;
      }

    } //for vars

    if( !foundVar )
      ofs << varName << " " << value << " " << error << std::endl;

  } // while ifs good

  ofs.close();

}
