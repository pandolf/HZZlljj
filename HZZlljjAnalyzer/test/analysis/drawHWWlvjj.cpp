#include <stdlib.h>
#include <iostream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"



void draw_vs_pt_plots( DrawBase* db, int nPtBins, Double_t* ptBins, const std::string& histoName, const std::string& axisName, const std::string& units="", const std::string& instanceName="Entries", bool log=false );


int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawHWWlvjj [(string) selectionType] [(string) leptType=\"ALL\"]" << std::endl;
    exit(23);
  }

  std::string selType(argv[1]);

 // std::string zJets_dataset = "ZJets_alpgen";

  std::string leptType = "ALL";
  if( argc==3 ) {
    std::string leptType_str(argv[2]);
    leptType = leptType_str;
  }

  DrawBase* db = new DrawBase("HWWlvjj");
  db->set_pdf_aussi((bool)false);

  std::string outputdir_str = "HWWlvjjPlots_"+selType;
  //if( leptType!="ALL" ) outputdir_str += "_" + leptType;
  outputdir_str += "_" + leptType;
  db->set_outputdir(outputdir_str);

   TFile* mcSignal_WW200 = TFile::Open("HWWlvjj_WW200_loose_ALL.root");
   db->add_mcFile(mcSignal_WW200, "WW200", "HWW (200)", 2);
 //  TFile* mcSignal_WW300 = TFile::Open("HWWlvjj_WW300_presel_ALL.root");
 // db->add_mcFile(mcSignal_WW300, "WW300", "HWW (300)", 46);
 //  TFile* mcSignal_WW400 = TFile::Open("HWWlvjj_WW400_opt400_ALL.root");
 // db->add_mcFile(mcSignal_WW400, "WW400", "HWW (400)", 93);
  //  TFile* mcSignal_WW500 = TFile::Open("HWWlvjj_WW500_loose_ALL.root");
  //  db->add_mcFile(mcSignal_WW500, "WW500", "HWW (500)", 98);

 // TFile* mcBkg_WJets = TFile::Open("HWWlvjj_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1_2_helicity_ALL.root");
 // db->add_mcFile(mcBkg_WJets, "WJets", "W + Jets", 38);
 // TFile* mcBkg_VV = TFile::Open("HWWlvjj_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_2_helicity_ALL.root");
 // db->add_mcFile(mcBkg_VV, "VV", " Diboson ", 65);
 // TFile* mcBkg_TT = TFile::Open("HWWlvjj_TT_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_3_helicity_ALL.root");
 // db->add_mcFile(mcBkg_TT, "tt", " tt ", 30);
  //TFile* mcBkg_DY = TFile::Open("HWWlvjj_DY_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1_2_helicity_ALL.root");
//  db->add_mcFile(mcBkg_DY, "DY", " DY ", 39);
 // TFile* mcBkg_T = TFile::Open("HWWlvjj_TToBLNu_TuneZ2_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1_2_helicity_ALL.root");
  //db->add_mcFile(mcBkg_T, "top", " top ", 50);


 // db->set_lumiNormalization( 1000. ); //1 fb-1 li somma e norm alla lumi
  //db->set_noStack( false);//(bool)true );
  db->set_shapeNormalization();


 // db->drawHisto( "nJets_presel", "Number of jets after preselection","", "Events",false);
 // db->drawHisto( "EtaOtherJets", "Eta Leading Jet not from H","", "Events",false);

 // db->drawHisto( "mtW_JustPresel", "W Boson Transverse Mass", "Gev/c^{2}", "Events", false);
 // db->set_rebin( 10 );
 // db->drawHisto( "mWW_kinfit", "H mass after all cuts", "Gev/c^{2}", "Events", false);
 // db->drawHisto( "helicityLD_kinfit", "Helicity Angles Likelihood Discriminant", "", "Events", false);
//  db->drawHisto( "QGLikelihoodProd", "Quark-Gluon Likelihood", "", "Events", false);
//  db->drawHisto( "QGLikelihoodJet1", "Jet1 Quark-Gluon Likelihood", "", "Events", false);
//  db->drawHisto( "QGLikelihoodJet2", "Jet2 Quark-Gluon Likelihood", "", "Events", false);

// db->drawHisto( "btag", "B-Tag Value for all jet after preselection", "", "Events", false);

//  db->drawHisto( "ptJet1_JustPresel", "Leading Jet Transverse Momentum", "Gev/c", "Events", false);
//  db->drawHisto( "ptJet2_JustPresel", "Sublead. Jet Transverse Momentum", "Gev/c", "Events", false);
//  db->drawHisto( "ptLept1_JustPresel", "Leading Lepton Transverse Momentum", "Gev/c", "Events", false);
//  db->drawHisto( "ptLept2_JustPresel", "Sublead. Lepton Transv. Momentum", "Gev/c", "Events", false);
//  db->drawHisto( "mWjj_JustPresel", "Dijet Invariant Mass", "Gev/c^{2}", "Events", false);
//  db->drawHisto( "deltaRll_JustPresel", "Lepton-Lepton #Delta R", "", "Events", false);
//  db->drawHisto( "deltaRjj_JustPresel", "Jet-Jet #Delta R", "", "Events", false);
//  db->drawHisto( "energyMet", "Missing Energy", "GeV", "Events", false);
   
   //db->set_rebin(2);
   //db->set_xAxisMin(400.);
   //db->set_yAxisMax(500);
   //db->drawHisto( "mWW_kinfit", "WW Inv. Mass After all cuts", "GeV/c^{2}", "Events", false, 2);

 // db->drawHisto( "ptJet1", "Pt Leading Jet", "GeV/c2", "Events", false);
 // db->drawHisto( "ptJet2", "Pt Leading Jet", "GeV/c2", "Events", false);
              // Resolution before and after Fit Jet1
/*  std::vector< HistoAndName > ptResoJet1_withFit;
  HistoAndName ptResoJet1_beforeKin;
  ptResoJet1_beforeKin.histoName = "ptResoJet1_beforeKin";
  ptResoJet1_beforeKin.legendName = "Before KinFit";
  ptResoJet1_withFit.push_back( ptResoJet1_beforeKin );
  HistoAndName ptResoJet1_afterKin;
  ptResoJet1_afterKin.histoName = "ptResoJet1_afterKin";
  ptResoJet1_afterKin.legendName = "After KinFit";
  ptResoJet1_withFit.push_back( ptResoJet1_afterKin ); 
  //char legendResoJet1Title[200];
  //db->set_legendTitle( legendResoJet1Title );
  db->compareDifferentHistos( ptResoJet1_withFit, "ptResoJet1_withFit", "Pt Resolution Jet1",  "");
              // Resolution before and after Fit Jet2
  std::vector< HistoAndName > ptResoJet2_withFit;
  HistoAndName ptResoJet2_beforeKin;
  ptResoJet2_beforeKin.histoName = "ptResoJet2_beforeKin";
  ptResoJet2_beforeKin.legendName = "Before KinFit";
  ptResoJet2_withFit.push_back( ptResoJet2_beforeKin );
  HistoAndName ptResoJet2_afterKin;
  ptResoJet2_afterKin.histoName = "ptResoJet2_afterKin";
  ptResoJet2_afterKin.legendName = "After KinFit";
  ptResoJet2_withFit.push_back( ptResoJet2_afterKin );
  db->compareDifferentHistos( ptResoJet2_withFit, "ptResoJet2_withFit", "Pt Resolution Jet2",  "");
              // Resolution before and after Fit ptW
  std::vector< HistoAndName > ptWReso_withFit;
  HistoAndName ptWReso_beforeKin;
  ptWReso_beforeKin.histoName = "ptWreso_beforeKin";
  ptWReso_beforeKin.legendName = "Before KinFit";
  ptWReso_withFit.push_back( ptWReso_beforeKin );
  HistoAndName ptWReso_afterKin;
  ptWReso_afterKin.histoName = "ptWreso_afterKin";
  ptWReso_afterKin.legendName = "After KinFit";
  ptWReso_withFit.push_back( ptWReso_afterKin );
  db->compareDifferentHistos( ptWReso_withFit, "ptWReso_withFit", "Pt Resolution W(jj)",  "");
              // Resolution before and after Fit M(H)
  std::vector< HistoAndName > mHReso_withFit;
  HistoAndName mHReso_beforeKin;
  mHReso_beforeKin.histoName = "mWW_nokinfit";
  mHReso_beforeKin.legendName = "Before KinFit";
  mHReso_withFit.push_back( mHReso_beforeKin );
  HistoAndName mHReso_afterKin;
  mHReso_afterKin.histoName = "mWW_kinfit";
  mHReso_afterKin.legendName = "After KinFit";
  mHReso_withFit.push_back( mHReso_afterKin );
  db->compareDifferentHistos( mHReso_withFit, "mHReso_withFit", "m(H) Resolution",  "");
             // M(H) con le due soluzioni
  std::vector< HistoAndName > mHrightSol;
  HistoAndName MH_rightsol;
  MH_rightsol.histoName = "mHrightSol";
  MH_rightsol.legendName = "M(H) right neu";
  mHrightSol.push_back( MH_rightsol );
  HistoAndName MH_wrongsol;
  MH_wrongsol.histoName = "mHwrongSol";
  MH_wrongsol.legendName = "M(H) wrong neu";
  mHrightSol.push_back( MH_wrongsol );
  db->compareDifferentHistos( mHrightSol, "mHiggs_2Sol", "m(H)",  "");
              // Eta(H) con le due soluzioni
  std::vector< HistoAndName > etaHrightSol;
  HistoAndName etaH_rightsol;
  etaH_rightsol.histoName = "etaHrightSol";
  etaH_rightsol.legendName = "#eta(H) right neu";
  etaHrightSol.push_back( etaH_rightsol );
  HistoAndName etaH_wrongsol;
  etaH_wrongsol.histoName = "etaHwrongSol";
  etaH_wrongsol.legendName = "#eta(H) wrong neu";
  etaHrightSol.push_back( etaH_wrongsol );
  db->compareDifferentHistos( etaHrightSol, "etaHiggs_2Sol", "#eta(H)",  "");
             // Resolution After Pz
  std::vector< HistoAndName > ResoMet;
  HistoAndName Reso_MC;
  Reso_MC.histoName = "resoPzGetMCRight;
  Reso_MC.legendName = "Inversion MC";
  ResoMet.push_back( Reso_MC );
  HistoAndName Reso_MET;
  Reso_MET.histoName = "resoPzGetRight;
  Reso_MET.legendName = "Inversion MET";
  ResoMet.push_back( Reso_MET );
  db->compareDifferentHistos( ResoMet, "Reso_MC-MET", "Neu Pz Resolution",  "");
             // Pz Reso GetPZMC RW
  std::vector< HistoAndName > GetPzRWMC;
  HistoAndName GetPz_pzRMC;
  GetPz_pzRMC.histoName = "pzResoNeut_GetPzMC_Right";
  GetPz_pzRMC.legendName = "Right sol";
  GetPzRWMC.push_back( GetPz_pzRMC );
  HistoAndName GetPz_pzWMC;
  GetPz_pzWMC.histoName = "pzResoNeut_GetPzMC_Wrong";
  GetPz_pzWMC.legendName = "Wrong sol";
  GetPzRWMC.push_back( GetPz_pzWMC );
  db->compareDifferentHistos( GetPzRWMC, "Pz_ResolMC_Get", "Neu MC Pz Resolution",  "");
             // Pz Reso GetPZ RW
  std::vector< HistoAndName > GetPzRW;
  HistoAndName GetPz_pzR;
  GetPz_pzR.histoName = "pzResoNeut_GetPz_Right";
  GetPz_pzR.legendName = "Right Sol";
  GetPzRW.push_back( GetPz_pzR );
  HistoAndName GetPz_pzW;
  GetPz_pzW.histoName = "pzResoNeut_GetPz_Wrong";
  GetPz_pzW.legendName = "Wrong Sol";
  GetPzRW.push_back( GetPz_pzW );
  HistoAndName GetPz_pzO;
  GetPz_pzO.histoName = "pzResoNeut_GetPz_One";
  GetPz_pzO.legendName = "One Sol";
  GetPzRW.push_back( GetPz_pzO );
  db->compareDifferentHistos( GetPzRW, "Pz_ResolRW_Get", "Neu Pz Resolution","","Entries",true);
          */   // Pz Reso RW
  std::vector< HistoAndName > FitRW;
  HistoAndName Fit_pzR;
  Fit_pzR.histoName = "Studio1";//"pzResoNeut_KinFit_Right";
  Fit_pzR.legendName = "No KinFit";//"Right Sol";
  FitRW.push_back( Fit_pzR );
  HistoAndName Fit_pzW;
  Fit_pzW.histoName = "Studio2";//"pzResoNeut_KinFit_Wrong";
  Fit_pzW.legendName = "KinFit Lept-Neu";//"Wrong Sol";
  FitRW.push_back( Fit_pzW);
  HistoAndName Fit_pzO;
  //Fit_pzO.histoName = "Studio3";//"pzResoNeut_KinFit_One";
  //Fit_pzO.legendName = "Global KinFit";//"One Sol";
  //FitRW.push_back( Fit_pzO);
  db->compareDifferentHistos( FitRW, "Higgsreso_fits","Pt Higgs Resolution","");//FitRW, "Pz_ResolRW_Fit", "Neu Pz Resolution","","Entries",true);
   /*         // Improvment on Pt Neu
  std::vector< HistoAndName > ImpPtNeu_Fit2;
  HistoAndName ptFit2;
  ptFit2.histoName = "ptResoNeut_KinFit";
  ptFit2.legendName = "KinFit";
  ImpPtNeu_Fit2.push_back( ptFit2 );
  HistoAndName ptGet;
  ptGet.histoName = "ptResoNeut_GetPz";
  ptGet.legendName = "No KinFit ";
  ImpPtNeu_Fit2.push_back( ptGet);
  db->compareDifferentHistos( ImpPtNeu_Fit2, "ImpPtNeu_Fit", "Pt Resolution",  "");
             // Improvment on Pz Neu
  std::vector< HistoAndName > ImpPzNeu_Fit2;
  HistoAndName pzFit2;
  pzFit2.histoName = "pzResoNeut_KinFit";
  pzFit2.legendName = "KinFit";
  ImpPzNeu_Fit2.push_back( pzFit2 );
  HistoAndName pzGet;
  pzGet.histoName = "pzResoNeut_GetPz";
  pzGet.legendName = "No KinFit ";
  ImpPzNeu_Fit2.push_back( pzGet);
  db->compareDifferentHistos( ImpPzNeu_Fit2, "ImpPzNeu_Fit2", "Pz Resolution",  "");
             // Improvment on Pt W
  std::vector< HistoAndName > ImpPtW_Fit2;
  HistoAndName ptFit2W;
  ptFit2W.histoName = "ptResoW_KinFit";
  ptFit2W.legendName = "KinFit";
  ImpPtW_Fit2.push_back( ptFit2W );
  HistoAndName ptGetW;
  ptGetW.histoName = "resoPtWll";
  ptGetW.legendName = "No KinFit ";
  ImpPtW_Fit2.push_back( ptGetW );
  db->compareDifferentHistos( ImpPtW_Fit2, "ImpPtW_Fit2", "Pt Resolution",  "");
             // Improvment on Pz Neu
  std::vector< HistoAndName > ImpPzW_Fit2;
  HistoAndName pzFit2W;
  pzFit2W.histoName = "pzResoW_KinFit";
  pzFit2W.legendName = "KinFit";
  ImpPzW_Fit2.push_back( pzFit2W );
  HistoAndName pzGetW;
  pzGetW.histoName = "resoPzWll";
  pzGetW.legendName = "No KinFit ";
  ImpPzW_Fit2.push_back( pzGetW );
  db->compareDifferentHistos( ImpPzW_Fit2, "ImpPzW_Fit2", "Pz Resolution",  "");
            // Chi-Square Fit/fitJet
  std::vector< HistoAndName > chisquare_FitJets;
  HistoAndName chisquare_Fit;
  chisquare_Fit.histoName = "kinfit2_chiSquareProb";
  chisquare_Fit.legendName = "Fit Lept-Neu";
  chisquare_FitJets.push_back( chisquare_Fit );
  HistoAndName chisquare_Jets;
  chisquare_Jets.histoName = "kinfit_chiSquareProb";
  chisquare_Jets.legendName = "Fit Jets";
  chisquare_FitJets.push_back( chisquare_Jets );
  db->compareDifferentHistos( chisquare_FitJets, "ChiSquare_Fit", "Chi Square Prob.",  "");

  db->drawHisto( "pzResoNeut_KinFit", "Pz Resolution", "", "Events", false);

             // Helicity WR
  std::vector< HistoAndName > helicityRW;
  HistoAndName helicityR;
  helicityR.histoName = "helicityLDRight";
  helicityR.legendName = "Right Helicity";
  helicityRW.push_back( helicityR );
  HistoAndName helicityW;
  helicityW.histoName = "helicityLDWrong";
  helicityW.legendName = "Wrong Helicity";
  helicityRW.push_back( helicityW );
  db->compareDifferentHistos( helicityRW, "helicityRW", "Helicity",  "");

  db->drawHisto( "diffHelicity", "Helicity Right-Helicity Wrong", "", "Events", false);

             // Reso Pz Neu with Helicity 
  std::vector< HistoAndName > resoPzNeuHeli;
  HistoAndName resoPzNeuH;
  resoPzNeuH.histoName = "pzResoNeut_heli";
  resoPzNeuH.legendName = "Using Helicity";
  resoPzNeuHeli.push_back( resoPzNeuH );
  HistoAndName resoPzNeuGet;
  resoPzNeuGet.histoName = "pzResoNeut_GetPz";
  resoPzNeuGet.legendName = "Using Invers.";
  resoPzNeuHeli.push_back( resoPzNeuGet );
  db->compareDifferentHistos( resoPzNeuHeli, "ResoPzNeu_helicity", "Neu Pz resolution",  "");
             // Reso Pz Neu with Helicity 
  std::vector< HistoAndName > resomHHeli;
  HistoAndName resomHH;
  resomHH.histoName = "ResomH_heli";
  resomHH.legendName = "Using Helicity";
  resomHHeli.push_back( resomHH );
  HistoAndName resomHGet;
  resomHGet.histoName = "ResomWW_GetPz";
  resomHGet.legendName = "Using Invers.";
  resomHHeli.push_back( resomHGet );
  db->compareDifferentHistos( resomHHeli, "ResomH_helicity", "M(H) resolution",  "");
            // Reso Pz Neu with Helicity 
  std::vector< HistoAndName > resoPzWHeli;
  HistoAndName resoPzWH;
  resoPzWH.histoName = "Studio3";//"ResoPzW_heli";
  resoPzWH.legendName = "Using Helicity";
  resoPzWHeli.push_back( resoPzWH );
  HistoAndName resoPzWGet;
  resoPzWGet.histoName = "Studio4";//"resoPzWll";
  resoPzWGet.legendName = "Using Invers.";
  resoPzWHeli.push_back( resoPzWGet );
  db->compareDifferentHistos( resoPzWHeli, "ResoPzW_helicity", "W Pz Resolution",  "");
*/
  delete db;
 
  db = 0;

  return 0;

}  


void draw_vs_pt_plots( DrawBase* db, int nPtBins, Double_t* ptBins, const std::string& histoName, const std::string& axisName, const std::string& units, const std::string& instanceName, bool log ) {

  for( unsigned iPtBin=0; iPtBin<nPtBins-1; ++iPtBin ) {

    Double_t ptMin = ptBins[iPtBin];
    Double_t ptMax = ptBins[iPtBin+1];

    char histoName_thisBin[250];
    sprintf( histoName_thisBin, "%s_pt_%.0lf_%.0lf", histoName.c_str(), ptMin, ptMax);

    std::string histoName_thisBin_str(histoName_thisBin);

    char legendTitle[200];
    sprintf( legendTitle, "%.0lf < pt < %.0lf GeV/c", ptMin, ptMax);
    std::string legendTitle_str(legendTitle);
    db->set_legendTitle(legendTitle_str);

    db->drawHisto( histoName_thisBin_str, axisName, units, instanceName, log );

  } //for bins

} //draw_vs_pt_plots

// To Have this Label without this code
//  TPaveText* cmsLabel = db->get_labelCMS();
//  TPaveText* sqrtLabel = db->get_labelSqrt();
//  cmsLabel->Draw("same"); 
//  sqrtLabel->Draw("same");
