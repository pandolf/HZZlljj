// @(#)root/tmva $Id: TMVAClassification.C 31458 2009-11-30 13:58:20Z stelzer $
/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l TMVAClassification.C\(\"Fisher,Likelihood\"\)                       *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set is used.                                     *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 **********************************************************************************/

#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>
#include <fstream>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Tools.h"
#endif



bool FIX_MZJJ = true;
bool PTZ = true;

// read input data file with ascii format (otherwise ROOT) ?
Bool_t ReadDataFromAsciiIFormat = kFALSE;
   
void TMVAClassification_helicityLD( int mass=400, int nbtags=-1 )
{
   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the 
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //
   // if you like to use a method via the plugin mechanism, we recommend using
   // 
   // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
   // (an example is given for using the BDT as plugin (see below),
   // but of course the real application is when you write your own
   // method based)

   // this loads the library
   TMVA::Tools::Instance();

   //---------------------------------------------------------------
   // default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   Use["Cuts"]            = 1;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // ---
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   // ---
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDERSkNN"]        = 0; // depreciated until further notice
   Use["PDEFoam"]         = 0;
   // --
   Use["KNN"]             = 0;
   // ---
   Use["HMatrix"]         = 0;
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0;
   Use["LD"]              = 0;
   // ---
   Use["FDA_GA"]          = 0;
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   // ---
   Use["MLP"]             = 0; // this is the recommended ANN
   Use["MLPBFGS"]         = 0; // recommended ANN with optional training method
   Use["CFMlpANN"]        = 0; // *** missing
   Use["TMlpANN"]         = 0; 
   // ---
   Use["SVM"]             = 0;
   // ---
   Use["BDT"]             = 0;
   Use["BDTD"]            = 0;
   Use["BDTG"]            = 0;
   Use["BDTB"]            = 0;
   // ---
   Use["RuleFit"]         = 0;
   // ---
   Use["Plugin"]          = 0;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

 //if (myMethodList != "") {
 //   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

 //   std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
 //   for (UInt_t i=0; i<mlist.size(); i++) {
 //      std::string regMethod(mlist[i]);

 //      if (Use.find(regMethod) == Use.end()) {
 //         std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
 //         for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
 //         std::cout << std::endl;
 //         return;
 //      }
 //      Use[regMethod] = 1;
 //   }
 //}


   //std::string mZjj_text = (FIX_MZJJ) ? "_FIXMZJJ" : "" ;
   //std::string ptZ_text = (PTZ) ? "_PTZ" : "" ;


   // Create a new root output file.
   char outfileName[200];
   //sprintf( outfileName, "TMVA_%d%s.root", mass, ptZ_text.c_str());
   if( nbtags>=0 )
     sprintf( outfileName, "TMVAHelFixQG_%d_btag%d.root", mass, nbtags );
   else
     sprintf( outfileName, "TMVAHelFixQG_%d.root", mass );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory will
   // then run the performance analysis for you.
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/ 
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in 
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile, 
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D" );

   // If you wish to modify default settings 
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
 //factory->AddVariable("ptLept1",  "Lead Lepton p_{T}", "GeV/c", 'F');
 //factory->AddVariable("absEtaLept1",  "Lead Lepton |eta|", "", 'F');
 //factory->AddVariable("ptLept2",  "Sublead Lepton p_{T}", "GeV/c", 'F');
///factory->AddVariable("absEtaLept2",  "Sublead Lepton |eta|", "", 'F');
 //factory->AddVariable("mZll",  "Dilepton Inv. Mass", "GeV/c^{2}", 'F');
 //factory->AddVariable("ptZll",  "Dilepton pt", "GeV", 'F');
 //factory->AddVariable("deltaRll",  "#DeltaR (Lepton-Lepton)", "", 'F');
 //factory->AddVariable( "ptJet1",  "Lead Jet p_{T}", "GeV/c", 'F');
///factory->AddVariable("absEtaJet1", "Lead Jet |eta|", "", 'F');
 //factory->AddVariable( "ptJet2",  "Sublead Jet p_{T}", "GeV/c", 'F');
///factory->AddVariable("absEtaJet2", "Sublead Jet |eta|", "", 'F');
 //factory->AddVariable( "ptJetRecoil",  "Recoil Jet p_{T}", "GeV/c", 'F');
   //factory->AddVariable( "QGLikelihoodJet1",  "Lead Jet QG LD", "", 'F');
   factory->AddVariable( "helicityLD",  "helicityLD", "", 'F');
   //factory->AddVariable( "QGLikelihoodJet1Jet2Recoil",  "QG LD 3Product", "", 'F');
   //factory->AddVariable( "QGLikelihoodJet2",  "Sublead Jet QG LD", "", 'F');
   //factory->AddVariable( "QGLikelihoodJetRecoil",  "Recoil Jet QG LD", "", 'F');
///factory->AddVariable("absEtaJetRecoil", "Recoil Jet |eta|", "", 'F');
///factory->AddVariable("deltaR_recoil_jet1", "#DeltaR (Recoil - Lead Jet)", "", 'F');
///factory->AddVariable("deltaR_recoil_Zjj", "#DeltaR (Recoil - Dijet)", "", 'F');
///factory->AddVariable("deltaR_recoil_Higgs", "#DeltaR (Recoil - Higgs)", "", 'F');
   //factory->AddVariable("mZjj",  "Dijet Inv. Mass", "GeV/c^{2}", 'F');
 //if( PTZ )
 //  factory->AddVariable("ptZjj",  "Dijet p_{T}", "GeV/c", 'F');
 //else
 //factory->AddVariable("deltaRjj",  "#DeltaR (Jet-Jet)", "", 'F');
 //factory->AddVariable("deltaRZZ",  "#DeltaR (Z-Z)", "", 'F');
///factory->AddVariable("deltaAbsEtaZZ",  "Z-Z Delta( |eta| )", "", 'F');
 //factory->AddVariable("absDeltaEtaZZ",  "Z-Z |Delta(eta)|", "", 'F');
///factory->AddVariable("absDeltaPhiZZ",  "Z-Z Delta #phi", "rad", 'F');
 //factory->AddVariable("ptZZ",  "Higgs Candidate p_{T}", "GeV/c", 'F');
 //factory->AddVariable("pfMet",  "PFMet", "GeV", 'F');
   factory->AddVariable("mZZ",  "ZZ Inv. Mass", "GeV/c^{2}", 'F');
 //factory->AddVariable("nBTags",  "Number of btagged jets", "", 'I');

   // You can add so-called "Spectator variables", which are not used in the MVA training, 
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the 
   // input variables, the response values of all trained MVAs, and the spectator variables
   //factory->AddSpectator( "mZjj",  "DiJet Inv. Mass", "GeV/c^{2}", 'F' );
   //factory->AddSpectator( "mZll",  "DiLepton Inv. Mass", "GeV/c^{2}", 'F' );
   //factory->AddSpectator( "mZZ",  "ZZ Inv. Mass", "GeV/c^{2}", 'F' );

   // read training and test data
    // load the signal and background event samples from ROOT trees


    char signalFileName[300];
    //sprintf(signalFileName, "TMVA_2ndLevelTreeW_SMHiggsToZZTo2L2Q_M-%d_7TeV-jhu-pythia6.root", mass);
    sprintf(signalFileName, "TMVA_2ndLevelTreeW_SMHiggsToZZTo2L2Q_M-%d_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1.root", mass);
    TFile* signalFile = TFile::Open(signalFileName);
    TTree *signalTree = (TTree*)signalFile->Get("reducedTree");
    float signalWeight_tmp;
    signalTree->SetBranchAddress("eventWeight", &signalWeight_tmp);
    signalTree->GetEntry(0);
    float signalWeight = signalWeight_tmp;

    std::vector<TFile*> backgroundFiles;
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z0Jets_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z1Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z1Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z1Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_4.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z2Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z2Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z2Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z2Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z3Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z3Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z3Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z3Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z4Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z4Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z4Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z4Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z5Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z5Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z5Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_Z5Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_ZBB0JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_ZBB1JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_ZBB2JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_ZBB3JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_ZCC0JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_ZCC1JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_ZCC2JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1.root") );
    backgroundFiles.push_back( TFile::Open("TMVA_2ndLevelTreeW_ZCC3JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3.root") );


    std::vector<TTree*> backgroundTrees;
    std::vector<float> backgroundWeights;
    for( unsigned iFile=0; iFile<backgroundFiles.size(); ++iFile ) {
      TTree* thisTree = (TTree*)backgroundFiles[iFile]->Get("reducedTree");
      backgroundTrees.push_back(thisTree);
      float weight;
      thisTree->SetBranchAddress("eventWeight", &weight);
      thisTree->GetEntry(0);
      backgroundWeights.push_back(weight);
    }


    // ====== register trees ====================================================
    //
    // the following method is the prefered one:
    // you can add an arbitrary number of signal or background trees

  //factory->AddSignalTree    ( signalTree,     1. );
  //for( unsigned iTree=0; iTree<backgroundTrees.size(); ++iTree ) {
  //  std::cout << backgroundWeights[iTree]  << std::endl;
  //  factory->AddBackgroundTree( backgroundTrees[iTree], 1. );
  //}

    factory->AddSignalTree    ( signalTree,     signalWeight     );
    //factory->AddSignalTree    ( signalTree,     1.     );
    for( unsigned iTree=0; iTree<backgroundTrees.size(); ++iTree ) {
      factory->AddBackgroundTree( backgroundTrees[iTree], backgroundWeights[iTree] );
    }

    // To give different trees for training and testing, do as follows:
    //    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
    //    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );

    // Use the following code instead of the above two or four lines to add signal and background 
    // training and test events "by hand"
    // NOTE that in this case one should not give expressions (such as "var1+var2") in the input 
    //      variable definition, but simply compute the expression before adding the event
    // 
    //    // --- begin ----------------------------------------------------------
    //    std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
    //    Float_t  treevars[4];
    //    for (Int_t ivar=0; ivar<4; ivar++) signal->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
    //    for (Int_t i=0; i<signal->GetEntries(); i++) {
    //       signal->GetEntry(i);
    //       for (Int_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
    //       // add training and test events; here: first half is training, second is testing
    //       // note that the weight can also be event-wise	
    //       if (i < signal->GetEntries()/2) factory->AddSignalTrainingEvent( vars, signalWeight ); 
    //       else                            factory->AddSignalTestEvent    ( vars, signalWeight ); 
    //    }
    //
    //    for (Int_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
    //    for (Int_t i=0; i<background->GetEntries(); i++) {
    //       background->GetEntry(i); 
    //       for (Int_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
    //       // add training and test events; here: first half is training, second is testing
    //       // note that the weight can also be event-wise	
    //       if (i < background->GetEntries()/2) factory->AddBackgroundTrainingEvent( vars, backgroundWeight ); 
    //       else                                factory->AddBackgroundTestEvent    ( vars, backgroundWeight ); 
    //    }
    //    // --- end ------------------------------------------------------------
    //
    // ====== end of register trees ==============================================
   
   
   // This would set individual event weights (the variables defined in the 
   // expression need to exist in the original TTree)
   //    for signal    : factory->SetSignalWeightExpression("weight1*weight2");
   //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");
   //factory->SetBackgroundWeightExpression("weight");

   // Apply additional cuts on the signal and background samples (can be different)
   //TCut mycuts = "mZZ>190. && absEtaLept1<2.1 && ptLept1>35.0324 && deltaRll<2.14868 && ptJet1>34.2011 && ptJet2>22.5958 && mZjj>56.9411 && mZjj<112.121"; 
   //TCut mycuts = "";
   //TCut mycuts = "absEtaLept1<2.1";
   TCut mycuts;

   mycuts += "ptLept1>40.";

   mycuts += "mZll>70.";
   mycuts += "mZll<110.";

   mycuts += "mZjj>75.";
   mycuts += "mZjj<105.";

   mycuts += "QGLikelihoodJet1Jet2>0.1";

   if( nbtags==0 )  mycuts += "nBTags==0";
   if( nbtags==1 )  mycuts += "nBTags==1";
   if( nbtags==2 )  mycuts += "nBTags==2";

   //mycuts += "leptType==1"; //keep only electrons to avoid BR bug


 //if( mass>=500 ) {
 //    mycuts  = "ptLept1>40.";
 //    mycuts += "ptZll>150.";
 //    mycuts += "ptJet1>70.";
 //    mycuts += "ptJet2>60.";
 //    mycuts += "deltaRjj<2.";
 //} else if( mass>=400 ) {
 //    mycuts  = "ptLept1>40.";
 //    mycuts += "ptZll>70.";
 //    mycuts += "ptJet1>70.";
 //    mycuts += "ptJet2>30.";
 //    mycuts += "deltaRjj<2.";
 //    mycuts += "mZZ<440.";
 //    mycuts += "mZZ>360.";
 //} else if(mass>=300) {       
 //    mycuts  = "ptLept1>40.";
 //    mycuts += "ptZll>60.";
 //    mycuts += "ptJet1>40.";
 //    mycuts += "ptJet2>30.";
 //    mycuts += "deltaRjj<2.8";
 //} else {
 //    mycuts  = "ptLept1>30.";
 //    mycuts += "ptJet1>40.";
 //    mycuts += "ptJet2>30.";
 // }


     //char mZZ_min[100];
     //char mZZ_max[100];
     //sprintf( mZZ_min, "mZZ>%d", (int)(mass*0.9));
     //sprintf( mZZ_max, "mZZ<%d", (int)(mass*1.1));
     //mycuts+= mZZ_min;
     //mycuts+= mZZ_max;



   // tell the factory to use all remaining events in the trees after training for testing:
 //factory->PrepareTrainingAndTestTree( mycuts, mycutb,
 //                                     "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // If no numbers of events are given, half of the events in the tree are used for training, and 
   // the other half for testing:
       //factory->PrepareTrainingAndTestTree( "", "SplitMode=random:!V" );  
       factory->PrepareTrainingAndTestTree( mycuts, "SplitMode=random:!V" );  
   // To also specify the number of testing events, use:
   //    factory->PrepareTrainingAndTestTree( mycut, 
   //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );  

   // ---- Book MVA methods
   //
   // please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"]) {

      std::string bookConditions;
      bookConditions = "H:!V:FitMethod=MC"; 

      bookConditions += ":VarProp[0]=FMax"; //helicityLD
      bookConditions += ":VarProp[1]=NotEnforced"; //ZZ mass

    //} else {

    //  bookConditions += ":VarProp[0]=FMax"; //pt lead lepton
    //  bookConditions += ":VarProp[1]=FMin"; //deltaR lept lept
    //  bookConditions += ":VarProp[2]=FMax"; // pt lead jet
    //  bookConditions += ":VarProp[3]=FMax"; // pt sublead jet
    //  if( !FIX_MZJJ ) {
    //    bookConditions += ":VarProp[4]=NotEnforced"; //jet-jet inv mass
    //    bookConditions += ":VarProp[5]=FMin"; // deltaR jet jet
    //  } else {
    //    bookConditions += ":VarProp[4]=FMin"; // deltaR jet jet
    //  }
    //  //bookConditions += ":VarProp[6]=NotEnforced"; // Z-Z inv mass

    //}
      //bookConditions += ":EffSel:SampleSize=100000000";
      if( nbtags==-1 )
      bookConditions += ":EffSel:SampleSize=1000000";
      else
      bookConditions += ":EffSel:SampleSize=10000000";

      factory->BookMethod( TMVA::Types::kCuts, "Cuts", 
                           bookConditions.c_str() );
                           //"!H:!V:FitMethod=MC:VarProp[0]=FMin:EffSel:SampleSize=200000" );
                           //"!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   }


   if (Use["CutsD"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsD", 
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsPCA", 
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"]) {

      std::string bookConditions;
      bookConditions = "H:!V:FitMethod=GA"; 
      bookConditions += ":VarProp[0]=FMax"; //pt lead lepton
      bookConditions += ":VarProp[1]=FMin"; //deltaR lept lept
      bookConditions += ":VarProp[2]=FMax"; // pt lead jet
      bookConditions += ":VarProp[3]=FMax"; // pt sublead jet
      bookConditions += ":VarProp[4]=NotEnforced"; //jet-jet inv mass
      bookConditions += ":VarProp[5]=FMin"; // deltaR jet jet
      //bookConditions += ":VarProp[6]=NotEnforced"; // Z-Z inv mass
      bookConditions += ":EffSel:Steps=30:Cycles=3:PopSize=1000:SC_steps=20:SC_rate=5:SC_factor=0.95";

      factory->BookMethod( TMVA::Types::kCuts, "CutsGA", bookConditions.c_str() );
                         //"H:!V:FitMethod=GA:CutRangeMin[0]=0:CutRangeMax[0]=3.2:VarProp[0]=FMin:EffSel:Steps=5:Cycles=3:PopSize=100:SC_steps=5:SC_rate=5:SC_factor=0.95" );


   }


   
   if (Use["CutsSA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
   
   // Likelihood
   if (Use["Likelihood"])
      factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood", 
                           "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" ); 

   // test the decorrelated likelihood
   if (Use["LikelihoodD"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD", 
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" ); 

   if (Use["LikelihoodPCA"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA", 
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 
 
   // test the new kernel density estimator
   if (Use["LikelihoodKDE"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE", 
                           "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

   // test the mixed splines and kernel density estimator (depending on which variable)
   if (Use["LikelihoodMIX"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX", 
                           "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

   // test the multi-dimensional probability density estimator
   // here are the options strings for the MinMax and RMS methods, respectively:
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );   
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );   
   if (Use["PDERS"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERS", 
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

   if (Use["PDERSkNN"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSkNN", 
                           "!H:!V:VolumeRangeMode=kNN:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

   if (Use["PDERSD"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSD", 
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

   if (Use["PDERSPCA"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSPCA", 
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

   // Multi-dimensional likelihood estimator using self-adapting phase-space binning
   if (Use["PDEFoam"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam", 
                           "H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0333:nActiveCells=500:nSampl=2000:nBin=5:CutNmin=T:Nmin=100:Kernel=None:Compress=T" );

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod( TMVA::Types::kKNN, "KNN", 
                           "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );
   // H-Matrix (chi2-squared) method
   if (Use["HMatrix"])
      factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V" ); 

   // Fisher discriminant   
   if (Use["Fisher"])
      factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=60:NsmoothMVAPdf=10" );

   // Fisher with Gauss-transformed input variables
   if (Use["FisherG"])
      factory->BookMethod( TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

   // Composite classifier: ensemble (tree) of boosted Fisher classifiers
   if (Use["BoostedFisher"])
      factory->BookMethod( TMVA::Types::kFisher, "BoostedFisher", "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2");

   // Linear discriminant (same as Fisher)
   if (Use["LD"])
      factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None" );

	// Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );
   
   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

   if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_SA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   if (Use["FDA_MT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   if (Use["FDA_MCMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MCMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
      factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=500:HiddenLayers=N+5:TestRate=10:EpochMonitoring" );

   if (Use["MLPBFGS"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=500:HiddenLayers=N+5:TestRate=10:TrainingMethod=BFGS:!EpochMonitoring" );


   // CF(Clermont-Ferrand)ANN
   if (Use["CFMlpANN"])
      factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  
  
   // Tmlp(Root)ANN
   if (Use["TMlpANN"])
      factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...
  
   // Support Vector Machine
   if (Use["SVM"])
      factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );
   
   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTG", 
                           "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT", 
                           "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
   
   if (Use["BDTB"]) // Bagging
      factory->BookMethod( TMVA::Types::kBDT, "BDTB", 
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD", 
                           "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );
   
   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );
   
   // For an example of the category classifier, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // As an example how to use the ROOT plugin mechanism, book BDT via
   // plugin mechanism
   if (Use["Plugin"]) {
         //
         // first the plugin has to be defined, which can happen either through the following line in the local or global .rootrc:
         //
         // # plugin handler          plugin name(regexp) class to be instanciated library        constructor format
         // Plugin.TMVA@@MethodBase:  ^BDT                TMVA::MethodBDT          TMVA.1         "MethodBDT(TString,TString,DataSet&,TString)"
         // 
         // or by telling the global plugin manager directly
      gPluginMgr->AddHandler("TMVA@@MethodBase", "BDT", "TMVA::MethodBDT", "TMVA.1", "MethodBDT(TString,TString,DataSet&,TString)");
      factory->BookMethod( TMVA::Types::kPlugins, "BDT",
                           "!H:!V:NTrees=400:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=50" );
   }

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

     std::cout << std::endl << std::endl;
     std::cout << "##############" << std::endl;
     std::cout << "#      1     #" << std::endl;
     std::cout << "##############" << std::endl;
   // Train MVAs using the set of training events
   factory->TrainAllMethods();

     std::cout << std::endl << std::endl;
     std::cout << "##############" << std::endl;
     std::cout << "#      2     #" << std::endl;
     std::cout << "##############" << std::endl;
   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

     std::cout << std::endl << std::endl;
     std::cout << "##############" << std::endl;
     std::cout << "#      3     #" << std::endl;
     std::cout << "##############" << std::endl;
   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();    

   // --------------------------------------------------------------

   if (Use["Cuts"] || Use["CutsGA"]) {

     std::string cutsType = (Use["Cuts"]) ? "" : "GA";
     std::string methodName = (Use["Cuts"]) ? "Cuts" : "CutsGA";
     TMVA::IMethod* method = (TMVA::IMethod*)factory->GetMethod(methodName.c_str());

     TMVA::MethodCuts* cuts = dynamic_cast<TMVA::MethodCuts*>(method);
     
     for( unsigned iEff=1; iEff<11; ++iEff ) {


       char cutsFileName[500];
       //sprintf( cutsFileName, "cuts%s_%d_Seff%d%s%s.txt", cutsType.c_str(), mass, 10*iEff, mZjj_text.c_str(), ptZ_text.c_str());
       if( nbtags>=0 )
         sprintf( cutsFileName, "CUTSHEL/cuts%sHelFixQG_btag%d_%d_Seff%d.txt", cutsType.c_str(), nbtags, mass, 10*iEff );
       else
         sprintf( cutsFileName, "CUTSHEL/cuts%sHelFixQG_%d_Seff%d.txt", cutsType.c_str(), mass, 10*iEff );

       ofstream ofs(cutsFileName);

       std::vector<Double_t> cutsMin, cutsMax;
       cuts->GetCuts((float)iEff*0.10, cutsMin, cutsMax);

     //if( iEff==0 )
     //  cuts->GetCuts(0.05, cutsMin, cutsMax);
     //else
     //  cuts->GetCuts((float)iEff/10., cutsMin, cutsMax);
     
       if( cutsMin.size() != cutsMax.size() ) {
         std::cout << "WARNING!!! cutsMin.size() != cutsMax.size() !!! Exiting!" << std::endl;
         exit( 177 );
       }
     
       for( unsigned iCut=0; iCut<cutsMin.size(); ++iCut) 
         ofs << factory->DefaultDataSetInfo().GetVariableInfo(iCut).GetInternalName() << " "  << cutsMin[iCut] << " " << cutsMax[iCut] << std::endl;
     
       ofs << "ptLept1 40. 1e+30" << std::endl;
       ofs << "QGLikelihoodJet1Jet2 0.1 1e+30" << std::endl;
       ofs << "mZll 70. 110." << std::endl;
       ofs << "mZjj 75. 105." << std::endl;

     }  // for eff

   } // if cuts GA
   
   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;      

   signalFile->Close();

   for(unsigned iBGFile=0; iBGFile<backgroundFiles.size(); ++iBGFile )
     backgroundFiles[iBGFile]->Close();

   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
