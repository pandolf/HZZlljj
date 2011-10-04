#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"

#include "fitTools.h"
#include "DrawBase.h"


int main()  {

  DrawBase* db = new DrawBase("alpha");


  TH1F::AddDirectory(kTRUE);

  TChain* chain = new TChain("tree_passedEvents");
  chain->Add("HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root/tree_passedEvents");
  chain->Add("HZZlljjRM_TT_TW_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root/tree_passedEvents");
  chain->Add("HZZlljjRM_VV_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root/tree_passedEvents");

  TTree* tree = chain->CopyTree("");

  const int nBins_pt = 10;
  Double_t ptBins[nBins_pt];

  fitTools::getBins_int( nBins_pt, ptBins, 183., 800.);
  
  TH1D* h1_mZZ_signal_0btag = new TH1D("mZZ_signal_0btag", "", nBins_pt-1, ptBins);
  TH1D* h1_mZZ_signal_1btag = new TH1D("mZZ_signal_1btag", "", nBins_pt-1, ptBins);
  TH1D* h1_mZZ_signal_2btag = new TH1D("mZZ_signal_2btag", "", nBins_pt-1, ptBins);

  TH1D* h1_mZZ_sidebands_0btag = new TH1D("mZZ_sidebands_0btag", "", nBins_pt-1, ptBins);
  TH1D* h1_mZZ_sidebands_1btag = new TH1D("mZZ_sidebands_1btag", "", nBins_pt-1, ptBins);
  TH1D* h1_mZZ_sidebands_2btag = new TH1D("mZZ_sidebands_2btag", "", nBins_pt-1, ptBins);

  h1_mZZ_signal_0btag->Sumw2();
  h1_mZZ_signal_1btag->Sumw2();
  h1_mZZ_signal_2btag->Sumw2();

  h1_mZZ_sidebands_0btag->Sumw2();
  h1_mZZ_sidebands_1btag->Sumw2();
  h1_mZZ_sidebands_2btag->Sumw2();

std::cout << tree->Project("mZZ_signal_0btag", "mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==0)") << std::endl;;
std::cout << tree->Project("mZZ_signal_1btag", "mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==1)") << std::endl;;
std::cout << tree->Project("mZZ_signal_2btag", "mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==2)") << std::endl;;

std::cout << tree->Project("mZZ_sidebands_0btag", "mZZ", "eventWeight*(((mZjj<75. && mZjj>60.) || (mZjj>105. && mZjj<140.)) && nBTags==0)") << std::endl;;
std::cout << tree->Project("mZZ_sidebands_1btag", "mZZ", "eventWeight*(((mZjj<75. && mZjj>60.) || (mZjj>105. && mZjj<140.)) && nBTags==1)") << std::endl;;
std::cout << tree->Project("mZZ_sidebands_2btag", "mZZ", "eventWeight*(((mZjj<75. && mZjj>60.) || (mZjj>105. && mZjj<140.)) && nBTags==2)") << std::endl;;

TFile* fileprova = TFile::Open("PROVA.root", "recreate");
fileprova->cd();
  h1_mZZ_signal_0btag->Write();
  h1_mZZ_signal_1btag->Write();
  h1_mZZ_signal_2btag->Write();

  h1_mZZ_sidebands_0btag->Write();
  h1_mZZ_sidebands_1btag->Write();
  h1_mZZ_sidebands_2btag->Write();
fileprova->Close();


  TH1D* h1_alpha_0btag = new TH1D(*h1_mZZ_signal_0btag);
  TH1D* h1_alpha_1btag = new TH1D(*h1_mZZ_signal_1btag);
  TH1D* h1_alpha_2btag = new TH1D(*h1_mZZ_signal_2btag);

  h1_alpha_0btag->Divide(h1_mZZ_sidebands_0btag);
  h1_alpha_1btag->Divide(h1_mZZ_sidebands_1btag);
  h1_alpha_2btag->Divide(h1_mZZ_sidebands_2btag);

  h1_alpha_0btag->SetName("alpha_0btag");
  h1_alpha_1btag->SetName("alpha_1btag");
  h1_alpha_2btag->SetName("alpha_2btag");

  h1_alpha_0btag->SetMarkerStyle(20);
  h1_alpha_1btag->SetMarkerStyle(21);
  h1_alpha_2btag->SetMarkerStyle(22);

  h1_alpha_0btag->SetMarkerSize(1.6);
  h1_alpha_1btag->SetMarkerSize(1.6);
  h1_alpha_2btag->SetMarkerSize(1.6);

  h1_alpha_0btag->SetMarkerColor(46);
  h1_alpha_1btag->SetMarkerColor(38);
  h1_alpha_2btag->SetMarkerColor(kGray+1);


  TH2D* h2_axes = new TH2D("axes", "", 10, 183., 800., 10, 0., 2.5);
  h2_axes->SetYTitle("#alpha Factor");
  h2_axes->SetXTitle("m_{lljj} [GeV]");

  TPaveText* labelCMS = db->get_labelCMS();
  TPaveText* labelSqrt = db->get_labelSqrt();

  TLine* line_one = new TLine( 183., 1., 800., 1. );

  TLegend* legend = new TLegend(0.2, 0.65, 0.5, 0.9);
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry(h1_alpha_0btag, "0-tag Category", "P");
  legend->AddEntry(h1_alpha_1btag, "1-tag Category", "P");
  legend->AddEntry(h1_alpha_2btag, "2-tag Category", "P");

  TCanvas* c1 = new TCanvas("c1", "", 600., 600.);
  c1->cd();
  h2_axes->Draw();
  line_one->Draw("same");
  h1_alpha_0btag->Draw("p same");
  h1_alpha_1btag->Draw("p same");
  h1_alpha_2btag->Draw("p same");
  legend->Draw("same");
  labelCMS->Draw("same");
  labelSqrt->Draw("same");
  gPad->RedrawAxis();

  c1->SaveAs("alpha.eps");

  return 0;

}
