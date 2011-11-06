#include "cstdlib"

#include "TGraphAsymmErrors.h"
#include "TString.h"

#include "CommonTools/DrawBase.h"


std::pair<TGraphAsymmErrors*,TGraphAsymmErrors*> get_expectedLimit( const std::string& data_dataset, const std::string& PUType );
float getMedian( TH1D* h1 );
float get_lowBound( TH1D* h1, float frac );
float get_upBound( TH1D* h1, float frac );



int main( int argc, char* argv[] ) {

  if( argc!=2 ) {
    std::cout << "USAGE: ./draw_golfcourse [(string)data_dataset]" << std::endl;
    exit(23);
  }
  

  std::string data_dataset(argv[1]);

  TString data_dataset_tstr(data_dataset);
  std::string PUType = "Run2011A";
  if( data_dataset=="HR11" )
    PUType = "HR11";
  if( data_dataset_tstr.BeginsWith("Run2011B") )
    PUType = "Run2011B";


  DrawBase* db = new DrawBase("Golfcourse");


  std::pair<TGraphAsymmErrors*,TGraphAsymmErrors*> graphs_expected = get_expectedLimit( data_dataset, PUType );
  TGraphAsymmErrors* graphExpected68 = graphs_expected.first;
  TGraphAsymmErrors* graphExpected95 = graphs_expected.second;
graphExpected68->Print();
graphExpected95->Print();


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  graphExpected68->Draw("APE");
  graphExpected95->Draw("APEsame");
  c1->SaveAs("Prova.eps");

  return 0;


}




std::pair<TGraphAsymmErrors*,TGraphAsymmErrors*> get_expectedLimit( const std::string& data_dataset, const std::string& PUType ) {


  TGraphAsymmErrors* graph68 = new TGraphAsymmErrors(0);
  TGraphAsymmErrors* graph95 = new TGraphAsymmErrors(0);



  
  std::vector<int> masses;
  masses.push_back(200);
  masses.push_back(250);
  masses.push_back(300);
  masses.push_back(350);
  masses.push_back(400);
  masses.push_back(440);
  masses.push_back(500);
  masses.push_back(540);
  masses.push_back(600);

  TH1F::AddDirectory(kTRUE);

  for( unsigned imass = 0; imass<masses.size(); ++imass ) {


    std::string crab_suffix = "_4"; //crab submission suffix

    char fileName[400];
    sprintf( fileName, "%d_%s%s/res/mergedToys.root", masses[imass], data_dataset.c_str(), crab_suffix.c_str() );
    TFile* file = TFile::Open(fileName);
    TTree* limitTree = (TTree*)file->Get("limit");

    TH1D *h1_limit = new TH1D("h1_limit","",2000,0.,15.);
    limitTree->Project( "h1_limit", "limit" );


    float expectedLimit = getMedian( h1_limit );

    float up95 = get_upBound( h1_limit, 0.95 );
    float up68 = get_upBound( h1_limit, 0.68 );

    float low95 = get_lowBound( h1_limit, 0.95 );
    float low68 = get_lowBound( h1_limit, 0.68 );

    graph68->SetPoint( imass, masses[imass], expectedLimit );
    graph95->SetPoint( imass, masses[imass], expectedLimit );

    graph68->SetPointError( imass, 0., 0., expectedLimit-low68, up68-expectedLimit );
    graph95->SetPointError( imass, 0., 0., expectedLimit-low95, up95-expectedLimit );

  } 

  TH1F::AddDirectory(kFALSE);

  std::pair<TGraphAsymmErrors*,TGraphAsymmErrors*> graphPair;
  graphPair.first = graph68;
  graphPair.second = graph95;

  return graphPair;

}



float getMedian( TH1D* h1 ) {

  return get_upBound( h1, 0. ); //believe me

}


float get_upBound( TH1D* h1, float frac ) {

  int nToys = h1->GetEntries();
  int nBins = h1->GetXaxis()->GetNbins();
  float hRange = h1->GetXaxis()->GetXmax() - h1->GetXaxis()->GetXmin();

  float threshFrac = (1.-frac)/2.; //double sided
 
  float sum = 0.;
  int iBin = nBins+1;
  while (sum<(float)nToys*threshFrac && iBin>0){
    iBin--;
    sum += h1->GetBinContent(iBin);
  }

  float upBound = (float)iBin*hRange/(float)nBins;
  return upBound;

}


float get_lowBound( TH1D* h1, float frac ) {

  int nToys = h1->GetEntries();
  int nBins = h1->GetXaxis()->GetNbins();
  float hRange = h1->GetXaxis()->GetXmax() - h1->GetXaxis()->GetXmin();

  float threshFrac = (1.-frac)/2.; //double sided
 
  float sum = 0.;
  int iBin = 0;
  while (sum<(float)nToys*threshFrac && iBin<nBins){
    iBin++;
    sum += h1->GetBinContent(iBin);
  }

  float lowBound = (float)iBin*hRange/(float)nBins;
  return lowBound;

}

