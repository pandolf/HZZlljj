#include "Riostream.h"
#include "string"
#include <vector>

#include "TROOT.h"
#include "TStyle.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

//float SBlowCut=75;
float SBlowCut=60;
float SBhighCut=130;
//float SBhighCut=105;

TStyle *tdrStyle ;
void setTDRStyle() ;
void sidebands(int myBTAG=0,int isElectron=1){

setTDRStyle() ;
  // gPad->SetFrameFillColor(kWhite);
  //  gROOT->ProcessLine(".x tdrstyle.C");
  // gROOT->ProcessLine("setTDRStyle()");
 const int nMass=10;
 float mass_points[nMass]={200.0,250.0,300.0,350.0,400.0,450.0,500.0, 600.0, 700.0, 800.0};

  vector<string> bkg_filenames;
  bkg_filenames.push_back("HZZlljjRM_TT_TW_TuneZ2_7TeV-pythia6-tauola_Spring11_v2_optLD_looseBTags_v2_ALL.root");	      
  bkg_filenames.push_back("HZZlljjRM_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11_v2_optLD_looseBTags_v2_ALL.root");
  //bkg_filenames.push_back("../HZZlljjRM_ZBB_alpgen_TuneZ2_Spring11_v2_optLD_looseBTags_v2_ALL.root");
  //bkg_filenames.push_back("../HZZlljjRM_ZJets_alpgen_TuneZ2_Spring11_v2_optLD_looseBTags_v2_ALL.root");
  bkg_filenames.push_back("HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root");
  string MCname="MADGRAPH";
  //  filenames.push_back(""); 
  const int nBkg=bkg_filenames.size();

  //load trees
  // TFile *fbkg[nBkg];
  // TTree *tbkg[nBkg];
  // for(int ibkg=0;ibkg<nBkg;ibkg++){
  //   fbkg[ibkg]=new TFile( (bkg_filenames.at(ibkg)).c_str(),"READ")  ;
  //   tbkg[ibkg]=(TTree*)   fbkg[ibkg]->Get("tree_passedEvent");
  // }

  TH1F *hmc_sr=new TH1F("mc_sigreg","MC (TOTAL) distribtuion of M_{ZZ} in signal region;M_{ZZ} [GeV];#evt / 10 GeV",65, 150.0, 800.0);
  TH1F *hmc_sb=new TH1F("mc_sbreg","MC (TOTAL) distribtuion of M_{ZZ} in sideband region;M_{ZZ} [GeV];#evt / 10 GeV",65, 150.0, 800.0);
  TH1F *halpha=new TH1F("alpha","Ratio SignalRegion / SidebandRegion from MC;M_{ZZ} [GeV];#evt / 10 GeV",65, 150.0, 800.0);
  TH1F *hdata_sr=new TH1F("data_sigreg","DATA (175 inv pb) distribtuion of M_{ZZ} in signal region;M_{ZZ} [GeV];#evt / 10 GeV",65, 150.0, 800.0);
  TH1F *hdata_sb=new TH1F("data_sbreg","DATA (175 inv pb) distribtuion of M_{ZZ} in sideband region;M_{ZZ} [GeV];#evt / 10 GeV",65, 150.0, 800.0);
  TH1F *hdata_pred=new TH1F("data_pred","DATA (175 inv pb) distribtuion of M_{ZZ} in signal region predicted from sidebands;M_{ZZ} [GeV];#evt / 10 GeV",65, 150.0, 800.0);

  hmc_sr->Sumw2();
  hmc_sb->Sumw2();
  //loop over bkg samples and sum them in the total one (weighting for lumi) 
  float mc_sr_integral=0.0;
  float mc_sb_integral=0.0;
  for(int ibkg=0;ibkg<nBkg;ibkg++){
    cout<<"Loading "<< (bkg_filenames.at(ibkg)).c_str()<<endl;
    TFile *fbkg=new TFile( (bkg_filenames.at(ibkg)).c_str(),"READ")  ;
    // fbkg->ls();  
    TTree *tbkg=(TTree*)   fbkg->Get("tree_passedEvents");
    cout<<"Tree name: "<<flush<<tbkg->GetName()<<endl;
    bool isSB=true;
    int nbtags, leptType;
    float mzz, ew, mjj;
    tbkg->SetBranchAddress("mZZ",&mzz);
    tbkg->SetBranchAddress("eventWeight",&ew);
    tbkg->SetBranchAddress("isSidebands",&isSB);
    tbkg->SetBranchAddress("nBTags",&nbtags);
    tbkg->SetBranchAddress("mZjj",&mjj);
    tbkg->SetBranchAddress("leptType",&leptType);
    cout<<"CP 1"<<endl;
    for(int i=0;i<tbkg->GetEntries();i++){
      tbkg->GetEntry(i);
      if(((mjj>SBlowCut && mjj<75)||(mjj<SBhighCut && mjj>105))&&nbtags==myBTAG && mzz>183){
	mc_sb_integral+=ew;
	hmc_sb->Fill(mzz,1000.0*ew);
      }else if(mjj>75 && mjj<105 && nbtags==myBTAG && mzz>183){
	mc_sr_integral+=ew;  
	hmc_sr->Fill(mzz,1000.0*ew);
      }
      /*
      if(mzz<160.0||nbtags!=myBTAG)continue;
      if(mzz>160.0&&mzz<800.0){
	if((mjj > SBlowCut && mjj <75) || (mjj < SBhighCut && mjj > 105))mc_sb_integral+=ew	;
	else if(mjj>75 && mjj <105) mc_sr_integral+=ew	;
      }

      if(((mjj > SBlowCut && mjj <75) || (mjj < SBhighCut && mjj > 105))){
	hmc_sb->Fill(mzz,1000.0*ew);
      }
      else if(mjj>75 && mjj <105){
	hmc_sr->Fill(mzz,1000.0*ew);
      }
      */
    }//end loop over tree entries
    delete tbkg;
    delete fbkg;
  }//end loop on bkg files

  halpha->Sumw2();
  halpha->Divide(hmc_sr,hmc_sb);

  string rootFileName[3]={"alpha_0btag.root","alpha_1btag.root","alpha_2btag.root"};

  //$%^^%&@$%^@$^%@$#%   check this out... don't forget!
  // ------------------------  smoothing out alpha  -------------------------------------

  double BinContent=0;
  double SmoothingThreshold=3.0;
  for(int iBin=1; iBin<halpha->GetNbinsX(); iBin++){
    if(halpha->GetBinContent(iBin)>SmoothingThreshold){
      if(iBin!=halpha->GetNbinsX()){
	BinContent = halpha->GetBinContent(iBin);
	halpha->SetBinContent(iBin,(halpha->GetBinContent(iBin+1)+halpha->GetBinContent(iBin-1))/2);     
	cout << "WARNING: anomalous bin contents ... " << endl;
	cout << "Changing bin " << iBin << " from " << BinContent << " to " << halpha->GetBinContent(iBin) << endl;
      }else if(iBin==halpha->GetNbinsX()){
	BinContent = halpha->GetBinContent(iBin);
	halpha->SetBinContent(iBin,(halpha->GetBinContent(iBin-1)+halpha->GetBinContent(iBin))/2);
	cout << "WARNING: anomalous bin contents ... " << endl;
	cout << "Changing last bin from " << BinContent << " to " << halpha->GetBinContent(iBin) << endl;
      }else if(iBin==1){
	BinContent = halpha->GetBinContent(iBin);
	halpha->SetBinContent(iBin,(halpha->GetBinContent(iBin+1)+halpha->GetBinContent(iBin))/2);
	cout << "WARNING: anomalous bin contents ... " << endl;
	cout << "Changing first bin from " << BinContent << " to " << halpha->GetBinContent(iBin) << endl;
      }
    }
  }
 
  // --------------------------------------------------------------------------------------
  TFile *f_alpha = new TFile(rootFileName[myBTAG].c_str(),"UPDATE");

  string temp = "alpha_"+MCname;
  halpha->Write(temp.c_str());
  
  TLegend *lmc=new TLegend(0.7,0.8,0.94,0.94);
  lmc->AddEntry(hmc_sb,"MC Sideband region","P");
  lmc->AddEntry(hmc_sr,"MC Signal Region","P");
  
  TCanvas *c1=new TCanvas("can_mc","CANMC1",800,800);
  c1->cd();
  hmc_sr->SetMarkerStyle(20);
  hmc_sr->SetMarkerSize(1.4);
  hmc_sr->SetMarkerColor(kBlack);
  hmc_sr->Scale(1/hmc_sr->Integral(1,60));
  hmc_sb->SetMarkerStyle(21);
  hmc_sb->SetMarkerSize(1.4);
  hmc_sb->SetMarkerColor(kRed);
  hmc_sb->Scale(1/hmc_sb->Integral(1,60));
  hmc_sb->Draw("PE0");
  
  temp = "0btag"+MCname+"SB";
  hmc_sb->Write(temp.c_str());
  temp = "0btag"+MCname+"SigReg";
  hmc_sr->Write(temp.c_str());
  
  hmc_sr->Draw("PE0sames");

  lmc->Draw();

  f_alpha->cd();
  hmc_sb->Scale(1/hmc_sb->Integral());
  temp = "SBdist_"+MCname;
  hmc_sb->Write(temp.c_str());
  hmc_sr->Scale(1/hmc_sr->Integral());
  temp = "Sigdist_"+MCname;
  hmc_sr->Write(temp.c_str());

  float ratio_alpha_integral=mc_sr_integral/mc_sb_integral;
  cout<<"Integral MC Signal [160, 800]: "<<mc_sr_integral<<"   Sideband "<<mc_sb_integral<<"   Ratio: "<< ratio_alpha_integral <<endl;
  TCanvas *c2=new TCanvas("can_alpha","CANALPHA1",800,800);
  c2->cd();
  halpha->Draw("histE0");

  //return;
  //take data
  TFile *fdata=new TFile( "HZZlljjRM_DATA_LP11_optLD_looseBTags_v2_ALL.root","READ")  ;
  TTree *tdata=(TTree*)fdata->Get("tree_passedEvents");
  bool isSB=true;
  int nbtags;
  float mzz, ew,mjj;
  int nDATAsb=0;
  int nDATAsr=0;
  tdata->SetBranchAddress("mZZ",&mzz);
  tdata->SetBranchAddress("mZjj",&mjj);
  tdata->SetBranchAddress("eventWeight",&ew);
  tdata->SetBranchAddress("isSidebands",&isSB);
  tdata->SetBranchAddress("nBTags",&nbtags);
  tdata->SetBranchAddress("leptType",&leptType);

  hdata_sr->Sumw2();
  hdata_sb->Sumw2();
  for(int i=0;i<tdata->GetEntries();i++){
    tdata->GetEntry(i);
    if(mzz<183.0||nbtags!=myBTAG)continue;
    if(((mjj > SBlowCut && mjj <75) || (mjj < SBhighCut && mjj > 105))&& isElectron == leptType){
      hdata_sb->Fill(mzz,ew);
      nDATAsb++;
    }
    else if(mjj>75 && mjj<105 && leptType==isElectron){
      hdata_sr->Fill(mzz,ew);
      nDATAsr++;
    }
  }//end loop over tree entries
  cout << "N EVTS IN DATA SB: " << nDATAsb << " -> EXPECTED EVTS IN SR: " << nDATAsb*mc_sr_integral/mc_sb_integral << 
    " | ACTUAL EVTS IN SR: " << nDATAsr << endl;
  hdata_pred->Sumw2();
  hdata_pred->Multiply(hdata_sb,halpha);
  TCanvas *c3=new TCanvas("can_data","CANDATA",800,800);
  c3->cd();
  hdata_sr->SetMarkerStyle(20);
  hdata_sr->SetMarkerSize(1.4);
  hdata_sr->SetMarkerColor(kBlack);
  hdata_sb->SetMarkerStyle(21);
  hdata_sb->SetMarkerSize(1.4);
  hdata_sb->SetMarkerColor(kRed);
  hdata_pred->SetMarkerStyle(22);
  hdata_pred->SetMarkerSize(1.4);
  hdata_pred->SetMarkerColor(kBlue);
  TLegend *ldata=new TLegend(0.7,0.8,0.94,0.94);
  ldata->AddEntry(hdata_sb,"DATA Sideband region","P");
  ldata->AddEntry(hdata_sr,"DATA Signal Region","P");
  ldata->AddEntry(hdata_pred,"DATA Signal Region from sidebands","P");
  hdata_sb->Draw("PE0");
  hdata_sr->Draw("PE0sames");
  hdata_pred->Draw("PE0sames");
  ldata->Draw();

  f_alpha->cd();
  hdata_sb->Scale(1/hdata_sb->Integral());
  temp="SBdist_data";
  hdata_sb->Write(temp.c_str());
  
  //calc integrals
  float eSBtot=0.0, eSRtot=0.0;
  for(int im=0;im<nMass;im++){
    float lowm=mass_points[im]*0.94;
    float him=mass_points[im]*1.10;
    float epred=0.0;
    float esb=0.0;
    float eobs=0.0;
  
    for(int i=0;i<tdata->GetEntries();i++){
      tdata->GetEntry(i);
      if(mzz<160.0||nbtags!=myBTAG)continue;
   

      if(im==0){
      if(mzz>183.0&&mzz<800.0 && isElectron==leptType){
	if(((mjj > SBlowCut && mjj <75) || (mjj < SBhighCut && mjj > 105)))	eSBtot++;
	else if(mjj>75 && mjj<105) eSRtot++;
      }
      }

      if(mzz<lowm||mzz>him)continue;
      if(((mjj > SBlowCut && mjj <75) || (mjj < SBhighCut && mjj > 105)) && isElectron==leptType){
	epred+=halpha->GetBinContent(halpha->FindBin(mzz));
	esb++;
	//	float alphaerr=halpha->GetBinError(halpha->FindBin(mzz));	
      }
      else if(mjj>75 && mjj<105 && isElectron==leptType){
	eobs++;
      }
    }//end loop over tree entries
    cout<<"M="<< mass_points[im]<<"   Integral in range ["<<lowm << " , "<<him<<"] : Predicted from Sidebands = "<<epred<<"+/-"<< 100.0/sqrt(esb)<<" %(stat.)   Observed in data= "<<eobs<<" +/- "<<100.0/sqrt(eobs) <<" % (stat.)"<<endl;
  }//end loop over masses
  cout<<"\n\nIntegrated over entire range [160, 800]: DATAinSB = "<<eSBtot<<"  ; ratio_alpha = "<<ratio_alpha_integral<<"  ---> Predicted: "<<eSBtot*ratio_alpha_integral<<"  vs Observed: "<<eSRtot<<endl;

  f_alpha->Close();


}//end main










void setTDRStyle() {
  tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  //  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  // tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  //tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.04, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.5);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.04, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}
