///////////////////////////////////////////////////////////////////////
//
//    FILE: PUWeight.h
//   CLASS: PUWeight
// AUTHORS: I. Gonzalez Caballero
//    DATE: 09/03/2011
//
///////////////////////////////////////////////////////////////////////
#include "PUWeight.h"

// ROOT Includes
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"

// C++ includes
#include <iostream>
using namespace std;


//Set DEBUGPUWEIGHT to 1 to get some debug information. Set it to 2 for more
//detail debug information.
#define DEBUGPUWEIGHT 2

#ifdef DEBUG
#define DEBUGPUWEIGHT 1
#endif


PUWeight::PUWeight(float luminosity, const char* year, const std::string& idealMCPU):
  fData(0),
  fMC(0),
  fWeight(0) {

  //Load Data histogram
  if (!LoadDataHistogram(luminosity, year))
    return;

  //No MC given. Take ideal MC
  IdealMCHistogram(idealMCPU);


  //Calculate Weight
  CalculateWeight();
}


PUWeight::PUWeight(const char* mcfolder, const char* mcproccess, 
		   float luminosity, const char* year):
  fData(0),
  fMC(0),
  fWeight(0) {

  //Load Data histogram
  if (!LoadDataHistogram(luminosity, year))
    return;
  
  //Load MC Histogram
  if (!LoadMCHistogram(mcfolder, mcproccess))
    return;



  //Calculate Weight
  CalculateWeight();
}




TH1F* PUWeight::LoadMCHistogram(const char* mcfolder, const char* mcproccess) {
#ifdef DEBUGPUWEIGHT
  cout << ">> Getting pileup for the MC " << mcproccess 
       << " inside " << mcfolder << "..." << endl;
#endif
  
  TString dsfile;
  dsfile.Form("http://www.hep.uniovi.es/jfernan/PUhistos/%s/%s.root", 
	      mcfolder, mcproccess);
#if (DEBUGPUWEIGHT > 1)
  cout << "   + Opening " << dsfile << endl;
#endif
  
  TFile* fds = TFile::Open(dsfile);
  if (!fds) {
    cerr << "ERROR [PUWeight]: Could not open file " << dsfile << "!"  << endl
	 << "                  Revise dataset name (" << mcproccess 
	 << ") or internet connection" << endl;
    return 0;
  }
  
  //Read dataset histogram...
#if (DEBUGPUWEIGHT > 1)
  cout << "   + Looking for histogram..." << endl;
#endif
  
  fMC = (TH1F*) fds->Get("htemp")->Clone("PU_MC");
  if (!fMC) {
    cerr << "ERROR [PUWeight]: Could not find histogram for dataset " << mcproccess << "!"
	 << endl;
    return 0;
  }
  fMC->SetDirectory(0);

  if (fMC->Integral() != 1) {
    cout << "NOTE [PUWeight]: MC histogram is not normalized to 1! Normalizing..."
	 << endl;
    fMC->Scale(1./fMC->Integral());
  }

  fds->Close();
  return fMC;
  
}

void PUWeight::SetMCHistogram(const TH1F* mcHisto) {

  fMC = (TH1F*)mcHisto->Clone();

  if (fMC->Integral() != 1) {
    cout << "NOTE [PUWeight]: MC histogram is not normalized to 1! Normalizing..."
	 << endl;
    fMC->Scale(1./fMC->Integral());
  }

  CalculateWeight();
  
}


void PUWeight::SetDataHistogram(const TH1F* dataHisto) {

#ifdef DEBUGPUWEIGHT
  cout << ">> Switching to new data histogram: " << dataHisto->GetName()
       << endl;
#endif
  fData = (TH1F*)dataHisto->Clone("PU_Data");

  if (fData->Integral() != 1) {
    cout << "NOTE [PUWeight]: Data histogram is not normalized to 1! Normalizing..."
	 << endl;
    fData->Scale(1./fData->Integral());
  }

  fData->SetDirectory(0);
  CalculateWeight();

}


TH1F* PUWeight::LoadDataHistogram(float luminosity, const char* year) {

#ifdef DEBUGPUWEIGHT
  cout << ">> Getting pileup for the " << luminosity << " pb-1 of data..." 
       << endl;
#endif
  
  TString dtfile;
  TFile* fdt = 0;
  if (luminosity > 0) {
    dtfile.Form("http://www.hep.uniovi.es/jfernan/PUhistos/Data%s/PUdata_%.1f.root", 
		year, luminosity);

  
#if (DEBUGPUWEIGHT > 1)
    cout << "   + Opening " << dtfile << endl;
#endif

    fdt = TFile::Open(dtfile);
    if (!fdt) {
      cerr << "NOTE [PUWeight]: Could not find file " << dtfile << "!"  << endl;
      cerr << "                 Trying default PU profile for data..." << endl;
    }
  }

  if (!fdt) {
    dtfile="http://www.hep.uniovi.es/jfernan/PUhistos/Data2011A/PUdata.root";

#if (DEBUGPUWEIGHT > 1)
    cout << "   + Opening " << dtfile << endl;
#endif

    fdt = TFile::Open(dtfile);
    if (!fdt) {
      cerr << "ERROR [PUWeight]: Could not find default profile in \"" 
	   << dtfile << "\"!"  << endl
	   << "                  Is your internet connection working?" << endl;
      return 0;
    }
  }
  
  //Read data histogram...
  fData = (TH1F*) fdt->Get("pileup")->Clone("PU_Data");
  if (!fData) {
    cerr << "ERROR [PUWeight]: Could not find histogram for data!" << endl;
    return 0;
  }
  
  fData->SetDirectory(0);
  
  if (fData->Integral() != 1) {
    cout << "NOTE [PUWeight]: Data histogram is not normalized to 1! Normalizing..."
	 << endl;
    fData->Scale(1./fData->Integral());
  }

  fdt->Close();

  return fData;
}


TH1F* PUWeight::CalculateWeight() {
  if (fData && fMC) {
    unsigned int nbins = fData->GetXaxis()->GetNbins();
    float xmin  = fData->GetXaxis()->GetXmin();
    float xmax  = fData->GetXaxis()->GetXmax();
    fWeight = new TH1F("PUWeight", "PU Weight", nbins, xmin, xmax);
    fWeight->SetDirectory(0);
    fWeight->Divide(fData, fMC);
  }
  else {
    cerr << "ERROR [PUWeight]: Something weird happened when trying to calculate the weights."
	 << endl 
	 << "                  I could not find the data and/or mc histograms!"
	 << endl;
  }

  return fWeight;
}

TH1F* PUWeight::IdealMCHistogram( const std::string& puType) {
  unsigned int nbins = 25;
  float xmin = -0.5;
  float xmax = 24.5;

//if (fData) {
//  nbins = fData->GetXaxis()->GetNbins();
//  xmin  = fData->GetXaxis()->GetXmin();
//  xmax  = fData->GetXaxis()->GetXmax();
//}


  fMC = new TH1F("PU_MC", "PU^{MC} Weight", nbins, xmin, xmax);
 
/*
  float   idealpu[]  = {0.0698146584, 0.0698146584, 0.0698146584, 
		     0.0698146584, 0.0698146584, 0.0698146584,
		     0.0698146584, 0.0698146584, 0.0698146584,
		     0.0698146584, 0.0698146584, 0.0630151648,
		     0.0526654164, 0.0402754482, 0.0292988928,
		     0.0194384503, 0.0122016783, 0.007207042,
		     0.004003637,  0.0020278322, 0.0010739954,
		     0.0004595759, 0.0002229748, 0.0001028162,
		     4.58337152809607E-05};
*/

  float   idealpu[25];

  if( puType=="Spring11_Flat10" ) {

      idealpu[0]   =   0.0698146584;
      idealpu[1]   =   0.0698146584;
      idealpu[2]   =   0.0698146584;
      idealpu[3]   =   0.0698146584;
      idealpu[4]   =   0.0698146584;
      idealpu[5]   =   0.0698146584;
      idealpu[6]   =   0.0698146584;
      idealpu[7]   =   0.0698146584;
      idealpu[8]   =   0.0698146584;
      idealpu[9]   =   0.0698146584;
      idealpu[10]   =   0.0698146584;
      idealpu[11]   =   0.0630151648;
      idealpu[12]   =   0.0526654164;
      idealpu[13]   =   0.0402754482;
      idealpu[14]   =   0.0292988928;
      idealpu[15]   =   0.0194384503;
      idealpu[16]   =   0.0122016783;
      idealpu[17]   =   0.007207042;
      idealpu[18]   =   0.004003637;
      idealpu[19]   =   0.0020278322;
      idealpu[20]   =   0.0010739954;
      idealpu[21]   =   0.0004595759;
      idealpu[22]   =   0.0002229748;
      idealpu[23]   =   0.0001028162;
      idealpu[24]   =   4.58337152809607E-05;

  } else if( puType=="Summer11_S4" ) {

      idealpu[0] =  0.14551;
      idealpu[1] =  0.0644453;
      idealpu[2] =  0.0696412;
      idealpu[3] =  0.0700311;
      idealpu[4] =  0.0694257;
      idealpu[5] =  0.0685655;
      idealpu[6] =  0.0670929;
      idealpu[7] =  0.0646049;
      idealpu[8] =  0.0609383;
      idealpu[9] =  0.0564597;
      idealpu[10] =  0.0508014;
      idealpu[11] =  0.0445226;
      idealpu[12] =  0.0378796;
      idealpu[13] =  0.0314746;
      idealpu[14] =  0.0254139;
      idealpu[15] =  0.0200091;
      idealpu[16] =  0.0154191;
      idealpu[17] =  0.0116242;
      idealpu[18] =  0.00846857;
      idealpu[19] =  0.00614328;
      idealpu[20] =  0.00426355;
      idealpu[21] =  0.00300632;
      idealpu[22] =  0.00203485;
      idealpu[23] =  0.00133045;
      idealpu[24] =  0.000893794;

  } else {

    std::cout << std::endl << "-----> WARNING!!!!!!!!!!" << std::endl;
    std::cout << "Unknown puType: '" << puType << "'. Using Spring11_Flat10 default." << std::endl;
      idealpu[0]  =  0.0698146584;
      idealpu[1]   =   0.0698146584;
      idealpu[2]   =   0.0698146584;
      idealpu[3]   =   0.0698146584;
      idealpu[4]   =   0.0698146584;
      idealpu[5]   =   0.0698146584;
      idealpu[6]   =   0.0698146584;
      idealpu[7]   =   0.0698146584;
      idealpu[8]   =   0.0698146584;
      idealpu[9]   =   0.0698146584;
      idealpu[10]   =   0.0698146584;
      idealpu[11]   =   0.0630151648;
      idealpu[12]   =   0.0526654164;
      idealpu[13]   =   0.0402754482;
      idealpu[14]   =   0.0292988928;
      idealpu[15]   =   0.0194384503;
      idealpu[16]   =   0.0122016783;
      idealpu[17]   =   0.007207042;
      idealpu[18]   =   0.004003637;
      idealpu[19]   =   0.0020278322;
      idealpu[20]   =   0.0010739954;
      idealpu[21]   =   0.0004595759;
      idealpu[22]   =   0.0002229748;
      idealpu[23]   =   0.0001028162;
      idealpu[24]   =   4.58337152809607E-05;

  }

  for (unsigned int i = 0; i < nbins; i++) {
    if (i < 25)
      fMC->Fill(i, idealpu[i]);
    else
      fMC->Fill(i, 0);
  }
  return fMC;
}
