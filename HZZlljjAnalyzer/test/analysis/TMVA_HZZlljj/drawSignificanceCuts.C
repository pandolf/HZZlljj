#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH2D.h"
#include "TChain.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include <iostream>
#include <fstream>
#include "cl95cms.c"




TH1F* getHistoPassingCuts( std::string histoName, TTree* tree, int nbtags, std::vector<std::string> names, std::vector<float> cutsMin, std::vector<float> cutsMax, float massMin, float massMax );


void drawSignificanceCuts(int mass=400, int nbtags=-1, bool KINcuts=false, bool fixQG=false, bool fixHel=false, std::string type="") {

  if( type!="" && type!="KIN" ) {
    std::cout << "Type: " << type << " unknown." << std::endl;
    exit(10101);
  }

  TStyle *simpleStyle = new TStyle("simpleStyle","");
  simpleStyle->SetCanvasColor(0);
  simpleStyle->SetPadColor(0);
  simpleStyle->SetFrameFillColor(0);
  simpleStyle->SetStatColor(0);
  simpleStyle->SetOptStat(0);
  simpleStyle->SetTitleFillColor(0);
  simpleStyle->SetCanvasBorderMode(0);
  simpleStyle->SetPadBorderMode(0);
  simpleStyle->SetFrameBorderMode(0);
  simpleStyle->cd();


  char signalChainName[300];
  sprintf( signalChainName, "TMVA_2ndLevelTreeW_SMHiggsToZZTo2L2Q_M-%d_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1.root", mass );
  TChain* signalChain = new TChain("reducedTree");
  signalChain->Add( signalChainName );

//TFile* signalFile = TFile::Open("TMVA_2ndLevelTreeW_JHUgen_HiggsSM400_2l2j_FASTSIM.root" );
//TTree* signalChain = (TTree*)signalFile->Get("reducedTree");

  std::cout << "-> Added signal tree 'TMVA_2ndLevelTreeW_JHUgen_HiggsSM400_2l2j_FASTSIM.root/reducedTree'" << std::endl;
  std::cout << "-> Signal tree has " << signalChain->GetEntries() << " entries." << std::endl;

  TFile* file_s = TFile::Open(signalChainName);
  TH1D* h1_nCounter = (TH1D*)file_s->Get("nCounter");
  float nTotal_s = h1_nCounter->GetBinContent(1);

  TChain* bgChain = new TChain("reducedTree");
  bgChain->Add("TMVA_2ndLevelTreeW_Z0Jets_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z1Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z1Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z1Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_4.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z2Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z2Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z2Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z2Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z3Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z3Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z3Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z3Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z4Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z4Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z4Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z4Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z5Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z5Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z5Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root");
  bgChain->Add("TMVA_2ndLevelTreeW_Z5Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root");
  bgChain->Add("TMVA_2ndLevelTreeW_ZBB0JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3.root");
  bgChain->Add("TMVA_2ndLevelTreeW_ZBB1JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root");
  bgChain->Add("TMVA_2ndLevelTreeW_ZBB2JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root");
  bgChain->Add("TMVA_2ndLevelTreeW_ZBB3JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root");
  bgChain->Add("TMVA_2ndLevelTreeW_ZCC0JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root");
  bgChain->Add("TMVA_2ndLevelTreeW_ZCC1JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root");
  bgChain->Add("TMVA_2ndLevelTreeW_ZCC2JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root");
  bgChain->Add("TMVA_2ndLevelTreeW_ZCC3JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3.root");



  std::cout << "-> Added all bg trees." << std::endl;
  std::cout << "-> Total bg tree has " << bgChain->GetEntries() << " entries." << std::endl;


  char significanceFileName[500];
  std::string prefix;
  prefix = "significanceFile";


  std::string KINtext = (KINcuts) ? "KIN" : "Hel";
  std::string fixQGtext = (fixQG) ? "FixQG" : "";
  std::string fixHeltext = (fixHel) ? "_FIXHel" : "";
  
  sprintf( significanceFileName, "%s%s%s%s%s_btag%d_%d.txt", prefix.c_str(), type.c_str(), KINtext.c_str(), fixQGtext.c_str(), fixHeltext.c_str(), nbtags, mass);
  ofstream ofs_sign(significanceFileName);
  ofs_sign << "S eff.\tSigma  \tUL   \ts(1fb-1)\tb(1fb-1)\ts(Number)\tb(Number)" << std::endl;

  TGraphErrors* graph = new TGraphErrors(0);
  TGraphErrors* graphUL = new TGraphErrors(0);
  TGraphErrors* graphUL_bg30 = new TGraphErrors(0);
  float signMax = 0.;
  float effMax = 0.;
  float ul_min = 999999999999999.;
  float effS_ul_min = 0.;


  for( unsigned iEff=1; iEff<10; ++iEff ) {

    char infileName[500];
    if( KINtext=="KIN" ) 
      sprintf( infileName, "CUTS%s/cuts%s%s%s_btag%d_%d_Seff%d.txt", KINtext.c_str(), type.c_str(), fixQGtext.c_str(), fixHeltext.c_str(), nbtags, mass, iEff*10);
    else
      sprintf( infileName, "CUTS%s/cuts%s%s%s%s_btag%d_%d_Seff%d.txt", KINtext.c_str(), type.c_str(), KINtext.c_str(), fixQGtext.c_str(), fixHeltext.c_str(), nbtags, mass, iEff*10);
    //sprintf( infileName, "CUTS%s/cuts%s_btag%d%s%s_%d_Seff%d.txt", KINtext.c_str(), type.c_str(), nbtags, fixQGtext.c_str(), fixHeltext.c_str(), mass, iEff*10);
    ifstream ifs(infileName);
    std::cout << "-> Opening Seff file: " << infileName << std::endl;
  
    std::vector<std::string> varNames;
    std::vector<float> cutsMin;
    std::vector<float> cutsMax;

    while( ifs.good() && !ifs.eof() ) {

      std::string varName;
      float cutMin, cutMax;

      ifs >> varName >> cutMin >> cutMax;

      varNames.push_back( varName );
      cutsMin.push_back( cutMin );
      cutsMax.push_back( cutMax );

    } //while file is good
  
    ifs.close();

    // eliminate last element (last line is read and is empty):
    varNames.pop_back();
    cutsMin.pop_back();
    cutsMax.pop_back();

    float deltaMass = 200.;
    float massMin = (float)mass-deltaMass;
    float massMax = (float)mass+deltaMass;

    TH1F* h1_signal = getHistoPassingCuts( "signal", signalChain, nbtags, varNames, cutsMin, cutsMax, massMin, massMax);
    TH1F* h1_bg     = getHistoPassingCuts( "bg",         bgChain, nbtags, varNames, cutsMin, cutsMax, massMin, massMax);

    h1_signal->SetFillColor( 46 );
    h1_bg->SetFillColor( 38 );

  //// compute s/sqrt(b) in mass+-10% region
  //float massmin = mass*0.9;
  //float massmax = mass*1.1;
  //int binmin = h1_signal->FindBin(massmin);
  //int binmax = h1_signal->FindBin(massmax);
  //float s = h1_signal->Integral(binmin, binmax);
  //float b = h1_bg->Integral(binmin, binmax);

    float s = h1_signal->Integral(0, h1_signal->GetNbinsX());
    float b = h1_bg->Integral(0, h1_bg->GetNbinsX());

  //  b /= 1.3;

    float significance = (b>=0.) ? s / sqrt(b+1.5) : 0.;
    //float significance = (b>=0.) ? s / sqrt(b) : 0.;
    if( s>0. && b==0. ) significance = 10.;


    float effS = (float)h1_signal->GetEntries()/nTotal_s;
    if( effS>effMax ) effMax = effS;

    float ul = 1000.*CLA( 1000., 0., effS, 0., b, 0. );
    float ul_bg30 = 0.;
    if( KINcuts ) ul_bg30 = 1000.*CLA( 1000., 0., effS, 0., b, 0.05*b );
    else ul_bg30 = 1000.*CLA( 1000., 0., effS, 0., b, 0.05*b );
    //float ul_bg30 = 1000.*CLA( 1000., 0., effS, 0., b, 0.3*b );
    //float ul = 0.;
    if( ul < ul_min ) {
      ul_min = ul;
      effS_ul_min = effS;
    }


    graph->SetPoint( iEff-1, 100.*effS, significance );
    graphUL->SetPoint( iEff-1, 100.*effS, ul );
    graphUL_bg30->SetPoint( iEff-1, 100.*effS, ul_bg30 );

    if( significance > signMax ) 
      signMax = significance;

    float yMax = h1_signal->GetMaximum() + h1_bg->GetMaximum();
    yMax*=1.5;

    THStack* stack = new THStack();
    stack->Add( h1_bg );
    stack->Add( h1_signal );

    TH2D* h2_axes = new TH2D("axes", "", 10, massMin, massMax, 10, 0., yMax);
    h2_axes->SetXTitle("ZZ Invariant Mass [GeV/c^{2}]");
    h2_axes->SetYTitle("Events / fb^{-1}");
    h2_axes->GetXaxis()->SetTitleOffset(1.1);
    h2_axes->GetYaxis()->SetTitleOffset(1.5);


    TLegend* legend = new TLegend(0.6, 0.75, 0.88, 0.88);
    legend->SetFillColor(0);
    legend->SetTextSize(0.035);
    legend->AddEntry( h1_signal, "Signal", "F");
    legend->AddEntry( h1_bg, "Background", "F");

    char canvasName[250];
    sprintf( canvasName, "INVMASSPLOTS/invMassPlot%s%s%s%s_nbtag%d_%d_Seff%d.eps", type.c_str(), KINtext.c_str(), fixQGtext.c_str(), fixHeltext.c_str(), nbtags, mass, iEff*10);

    TPaveText* label = new TPaveText( 0.15, 0.65, 0.45, 0.85, "brNDC");
    label->SetFillColor(0);
    label->SetTextSize(0.035);
    label->AddText("L = 1 fb^{-1}");
    char signalLabel[100];
    sprintf( signalLabel, "s = %.2f (%d%%)", s, (int)(((float)h1_signal->GetEntries()/nTotal_s)*100) );
    label->AddText( signalLabel );
    char bgLabel[100];
    sprintf( bgLabel, "b = %.2f", b);
    label->AddText( bgLabel );
    char signifLabel[100];
    sprintf( signifLabel, "s / #sqrt{b} = %.2f", significance);
    label->AddText( signifLabel );

    
    TCanvas* c1 = new TCanvas("c1", "c1", 600., 600.);
    c1->cd();
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);
    h2_axes->Draw();
    stack->Draw("histo same");
    legend->Draw("same");
    label->Draw("same");
    //gPad->RedrawAxis();
    c1->SaveAs(canvasName);

    delete c1;
    delete legend;
    delete h2_axes;
    delete stack;
    

    ofs_sign << effS << "\t" << significance << "\t" << ul << "\t" << s << "\t" << b << "\t" << h1_signal->GetEntries() << "\t" << h1_bg->GetEntries() << std::endl;

    delete h1_signal;
    delete h1_bg;

std::cout << "### " << iEff << std::endl;
  } // for iEff

  std::cout << "> > >   BEST UL: " << ul_min << std::endl;
  std::cout << "> > >   signal eff: " << effS_ul_min << std::endl;

  ofs_sign.close();

  graph->SetMarkerSize(2.);
  graph->SetMarkerStyle(29);
  graph->SetMarkerColor(kRed+3);

  graphUL->SetMarkerSize(2.);
  graphUL->SetMarkerStyle(29);
  graphUL->SetMarkerColor(kRed+3);

  graphUL_bg30->SetMarkerSize(2.);
  graphUL_bg30->SetMarkerStyle(20);
  graphUL_bg30->SetMarkerColor(kOrange+1);

  TH2D* h2_axes_gr = new TH2D("axes_gr", "", 10, 0., 1.3*effMax*100., 10, 0., 1.6*signMax ); 
  //TH2D* h2_axes_gr = new TH2D("axes_gr", "", 10, 0., 1., 10, 0., 5.);
  h2_axes_gr->SetYTitle("S / #sqrt{B} (1 fb^{-1})");
  h2_axes_gr->SetXTitle("Signal Efficiency [%]");
  h2_axes_gr->GetXaxis()->SetTitleOffset(1.1);
  h2_axes_gr->GetYaxis()->SetTitleOffset(1.5);

  TH2D* h2_axes_grUL = new TH2D("axes_grUL", "", 10, 0., 1.3*effMax*100., 10, 0., 1300. );
  h2_axes_grUL->SetYTitle("UL [fb]");
  h2_axes_grUL->SetXTitle("Signal Efficiency [%]");
  h2_axes_grUL->GetXaxis()->SetTitleOffset(1.1);
  h2_axes_grUL->GetYaxis()->SetTitleOffset(1.5);


  TCanvas* c_gr = new TCanvas("c_gr", "c_gr", 600., 600.);
  c_gr->SetLeftMargin(0.12);
  c_gr->SetBottomMargin(0.12);
  c_gr->cd();

  TPaveText* mass_label = new TPaveText(0.5, 0.7, 0.8, 0.8, "brNDC");
  mass_label->SetFillColor(0);
  mass_label->SetTextSize(0.035);
  char mass_label_text[100];
  sprintf( mass_label_text, "m_{H} = %d GeV/c^{2}", mass);
  mass_label->AddText(mass_label_text);
  
  h2_axes_gr->Draw();
  graph->Draw("P same");
  mass_label->Draw("same");

  char sign_vs_Seff_name[250];
  sprintf(sign_vs_Seff_name, "significance%s%s%s%s_vs_Seff_nbtag%d_%d.eps", type.c_str(), KINtext.c_str(), fixQGtext.c_str(), fixHeltext.c_str(), nbtags, mass);
  c_gr->SaveAs(sign_vs_Seff_name);

  c_gr->Clear();

  TLegend* legend = new TLegend(0.5, 0.66, 0.88, 0.88, mass_label_text);
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( graphUL, "No Error on BG", "P" );
  if( KINcuts ) legend->AddEntry( graphUL_bg30, "10% Error on BG", "P" );
  else legend->AddEntry( graphUL_bg30, "5% Error on BG", "P" );

  h2_axes_grUL->Draw();
  legend->Draw("same");
  graphUL_bg30->Draw("P same");
  graphUL->Draw("P same");

  char ul_vs_Seff_name[250];
  sprintf(ul_vs_Seff_name, "UL%s%s%s%s_vs_Seff_nbtag%d_%d.eps", type.c_str(), KINtext.c_str(), fixQGtext.c_str(), fixHeltext.c_str(), nbtags, mass);
  c_gr->SaveAs(ul_vs_Seff_name);

}



TH1F* getHistoPassingCuts( std::string histoName, TTree* tree, int nbtags, std::vector<std::string> names, std::vector<float> cutsMin, std::vector<float> cutsMax, float massMin, float massMax ) {


  //std::cout << "::getHistoPassingCuts:: Begin." << std::endl;

  TH1F* histo = new TH1F(histoName.c_str(), "", 50, massMin, massMax);
  histo->Sumw2();


  Float_t eventWeight;
  tree->SetBranchAddress( "eventWeight", &eventWeight );

  Float_t absEtaLept1;
  tree->SetBranchAddress( "absEtaLept1", &absEtaLept1 );

  std::vector<float> variables(names.size());
  int index_mZZ=-1;
  int index_mZjj=-1;
  int index_mZll=-1;
  int index_ptLept1=-1;


  for( unsigned i=0; i<names.size(); ++i ) {
    //std::cout << "::getHistoPassingCuts:: Setting Branch Address: '" << names[i] << "'" << std::endl;
    if( names[i]=="QGLikelihoodJet1_T_QGLikelihoodJet2" ) continue;
  //  tree->SetBranchAddress("QGLikelihoodJet1"
  //} else {
      tree->SetBranchAddress( names[i].c_str(), &(variables[i]) );
//std::cout << "set " << names[i] << std::endl;
  //}
    if( names[i]=="mZZ" )
      index_mZZ = i;
    if( names[i]=="mZll" )
      index_mZll = i;
    if( names[i]=="mZjj" )
      index_mZjj = i;
    if( names[i]=="ptLept1" )
      index_ptLept1 = i;
  }

  Int_t leptType;
  tree->SetBranchAddress( "leptType", &leptType );

  Int_t nBTags;
  tree->SetBranchAddress( "nBTags", &nBTags );

  Float_t mZZ;
  if( index_mZZ<0 )
    tree->SetBranchAddress( "mZZ", &mZZ );

  Float_t mZjj;
  if( index_mZjj<0 )
    tree->SetBranchAddress( "mZjj", &mZjj );

  Float_t mZll;
  if( index_mZll<0 )
    tree->SetBranchAddress( "mZll", &mZll );

  Float_t ptLept1;
  if( index_ptLept1<0 )
    tree->SetBranchAddress( "ptLept1", &ptLept1 );


  int nentries = tree->GetEntries();

  //std::cout << "::getHistoPassingCuts:: Begin Loop." << std::endl;

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry) {
 
    tree->GetEntry(iEntry);

    bool pass = true;

    for( unsigned iVar=0; iVar<variables.size() /*&& pass*/; ++iVar) {

 //   // preselection:
 //   if( absEtaLept1>2.1 ) pass=false;

 //   if( names[iVar]=="mZZ" ) {
 //     if( variables[iVar]<190. ) pass=false;
 //   }
 //   else if( names[iVar]=="ptLept1" ) {
 //     if( variables[iVar]<35.0324 ) pass=false;
 //   }
 //   else if( names[iVar]=="deltaRll" ) {
 //     if( variables[iVar]>2.14868 ) pass=false;
 //   }
 //   else if( names[iVar]=="ptJet1" ) {
 //     if( variables[iVar]<34.2011 ) pass=false;
 //   }
 //   else if( names[iVar]=="ptJet2" ) {
 //     if( variables[iVar]<22.5958 ) pass=false;
 //   }
 //   else if( names[iVar]=="mZjj" ) {
 //     if( variables[iVar]<56.9411 || variables[iVar]>112.121 ) pass=false;
 //   }


      // selection:
//std::cout << "requiring " << names[iVar] << ">" << cutsMin[iVar] << " && " << names[iVar] << "<" << cutsMax[iVar] << std::endl;
      if( variables[iVar]<cutsMin[iVar] || variables[iVar]>cutsMax[iVar] )
        pass = false;

    }

    //if( leptType != 1 ) pass = false;

    if( index_ptLept1<0 )
      if( ptLept1<40. ) pass=false;

    if( nBTags != nbtags ) pass = false;


    if( !pass ) continue;

    if( index_mZZ<0 )
      histo->Fill( mZZ, 1000.*eventWeight ); //1 fb-1
    else
      histo->Fill( variables[index_mZZ], 1000.*eventWeight ); //1 fb-1

  } // for entries

  //std::cout << "::getHistoPassingCuts:: End Loop." << std::endl;

  return histo;

} // getHistoPassingCuts
