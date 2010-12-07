#include <stdlib.h>
#include <iostream>
#include <string>
#include "DrawBase.h"
#include "fitTools.h"





int main(int argc, char* argv[]) {



  DrawBase* db = new DrawBase("QG");

  TFile* mcFile = TFile::Open("QG_JHUgen_HiggsSM400_2l2j_FASTSIM_QG.root");
  db->add_mcFile(mcFile, "JHUgen_HiggsSM400_2l2j_FASTSIM_QG", "HZZlljj (400)", kRed+1);

  db->set_outputdir("QG_Plots");

  db->set_shapeNormalization();
  db->set_rebin(2);

  std::vector< HistoAndName > ptJet;
  HistoAndName hn1;
  hn1.histoName = "ptJet_gluon";
  hn1.legendName = "Gluon Jets";
  ptJet.push_back( hn1);
  HistoAndName hn2;
  hn2.histoName = "ptJet_quark";
  hn2.legendName = "Quark Jets";
  ptJet.push_back( hn2);
  db->compareDifferentHistos( ptJet, "ptJet", "Jet Transverse Momentum", "GeV/c", "Jets");

  std::vector< HistoAndName > ptD;
  hn1.histoName = "ptD_gluon";
  hn1.legendName = "Gluon Jets";
  ptD.push_back( hn1);
  hn2.histoName = "ptD_quark";
  hn2.legendName = "Quark Jets";
  ptD.push_back( hn2);
  db->compareDifferentHistos( ptD, "ptD", "Sum of Constituent p_{T}^{2} / Jet p_{T}^{2}", "", "Jets");

  std::vector< HistoAndName > rmsCand;
  hn1.histoName = "rmsCand_gluon";
  hn1.legendName = "Gluon Jets";
  rmsCand.push_back( hn1);
  hn2.histoName = "rmsCand_quark";
  hn2.legendName = "Quark Jets";
  rmsCand.push_back( hn2);
  db->compareDifferentHistos( rmsCand, "rmsCand", "Constituent RMS about Jet Axis", "", "Jets");

  std::vector< HistoAndName > nCharged;
  hn1.histoName = "nCharged_gluon";
  hn1.legendName = "Gluon Jets";
  nCharged.push_back( hn1);
  hn2.histoName = "nCharged_quark";
  hn2.legendName = "Quark Jets";
  nCharged.push_back( hn2);
  db->compareDifferentHistos( nCharged, "nCharged", "Jet Charged Multiplicity", "", "Jets");

  std::vector< HistoAndName > nNeutral;
  hn1.histoName = "nNeutral_gluon";
  hn1.legendName = "Gluon Jets";
  nNeutral.push_back( hn1);
  hn2.histoName = "nNeutral_quark";
  hn2.legendName = "Quark Jets";
  nNeutral.push_back( hn2);
  db->compareDifferentHistos( nNeutral, "nNeutral", "Jet Neutral Multiplicity", "", "Jets");



  delete db;
  db = 0;

  return 0;

}  


