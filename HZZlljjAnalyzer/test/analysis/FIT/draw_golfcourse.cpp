#include "cstdlib"

#include "TGraphAsymmErrors.h"
#include "TString.h"

#include "CommonTools/DrawBase.h"


std::pair<TGraphAsymmErrors*,TGraphAsymmErrors*> get_expectedLimit( const std::string& data_dataset, const std::string& PUType );
float getMedian( TH1D* h1 );
float get_lowBound( TH1D* h1, float frac );
float get_upBound( TH1D* h1, float frac );
TPaveText* get_labelCMS( DrawBase* db );



int main( int argc, char* argv[] ) {

  if( argc!=2 ) {
    std::cout << "USAGE: ./draw_golfcourse [(string)data_dataset]" << std::endl;
    exit(23);
  }
  

  std::string data_dataset(argv[1]);

  TString data_dataset_tstr(data_dataset);
  std::string PUType = "Run2011A";
  float lumi = 2100.;
  if( data_dataset=="HR11" ) {
    PUType = "HR11";
    lumi = 4200.;
  } if( data_dataset_tstr.BeginsWith("Run2011B") ) {
    PUType = "Run2011B";
  }

  DrawBase* db = new DrawBase("Golfcourse");
  db->set_lumiNormalization(lumi);
  TFile* file_data = TFile::Open("HZZlljjRM_DATA_Run2011A_FULL_optLD_looseBTags_v2_ALL.root"); //any file is good, got to fix this
  db->add_dataFile( file_data, "dummydata" );



  std::pair<TGraphAsymmErrors*,TGraphAsymmErrors*> graphs_expected = get_expectedLimit( data_dataset, PUType );
  TGraphAsymmErrors* graphExpected68 = graphs_expected.first;
  TGraphAsymmErrors* graphExpected95 = graphs_expected.second;

  graphExpected68->SetFillColor(kGreen);
  graphExpected68->SetLineColor(kGreen+2);
  graphExpected68->SetLineWidth(2);
  graphExpected68->SetFillStyle(1001);//solid

  graphExpected95->SetFillColor(kYellow);
  graphExpected95->SetFillStyle(1001);//solid

  TH2D* axes = new TH2D("axes", "", 10, 160., 640., 10, 0., 10.);
  axes->SetXTitle("m_{H} [GeV]");
  axes->SetYTitle("#sigma_{95%} / #sigma_{SM}");

  //TPaveText* labelCMS = get_labelCMS(db);
  TPaveText* labelCMS = db->get_labelCMS();
  TPaveText* labelSqrt = db->get_labelSqrt();

  TLine* line_one = new TLine( 160., 1., 640., 1.);
  line_one->SetLineColor(kRed);
  line_one->SetLineWidth(2);

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  axes->Draw();
  labelCMS->Draw("same");
  labelSqrt->Draw("same");
  graphExpected95->Draw("3same");
  graphExpected68->Draw("3same");
  graphExpected68->Draw("LXsame");
  line_one->Draw("same");
  gPad->RedrawAxis();
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



TPaveText* get_labelCMS( DrawBase* db ) {

  float x1 = 0.145;
  float y1 = 0.953;
  float x2 = 0.6;
  float y2 = 0.975;

  
  TPaveText* cmslabel = new TPaveText( x1, y1, x2, y2, "brNDC" );
  cmslabel->SetFillColor(kWhite);
  cmslabel->SetTextSize(0.038);
  cmslabel->SetTextAlign(11);
  cmslabel->SetTextSize(0.038);
  cmslabel->SetTextFont(62);
  std::string leftText = "CMS Preliminary 2011";
  cmslabel->SetTextAlign(11); // align left
  std::string lumiText = db->get_lumiText();
  cmslabel->AddText(Form("%s, %s", leftText.c_str(), lumiText.c_str()));

  return cmslabel;

}
