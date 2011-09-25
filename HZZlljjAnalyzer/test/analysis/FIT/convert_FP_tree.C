void convert_FP_tree(){

  // TFile *dfile = new TFile("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/bonato/HtoZZto2L2J-crab/trees_FP_Summer11/HZZlljjRM_DATA_800pb_optLD_looseBTags_v2_ALL_FP.root");
  // TFile *dfile = new TFile("/afs/cern.ch/user/w/whitbeck/scratch0/HZZlljjRM_DATA_1fb_optLD_looseBTags_v2_ALL.root");//859
  TFile *dfile = new TFile("treeFromFrancesco/HZZlljjRM_DATA_LP11_optLD_looseBTags_v2_ALL.root","READ");//1000
 TTree *t=(TTree*)dfile->Get("tree_passedEvents");
 cout<<"Loaded tree "<<t->GetName()<<endl;
 int leptTypeIN, nBTagsIN;
 float mZZIN, mZjjIN,eventWeightIN;
 cout<<"addressing"<<endl;
 t->SetBranchAddress("mZZ", &mZZIN);
 t->SetBranchAddress("mZjj", &mZjjIN);
 t->SetBranchAddress("nBTags", &nBTagsIN);
 t->SetBranchAddress("leptType", &leptTypeIN);
 t->SetBranchAddress("eventWeight", &eventWeightIN);


 double leptTypeOUT, nBTagsOUT;
 double mZZOUT, mZjjOUT,eventWeightOUT;
 TFile *fout=new TFile("./convertedTree_LP_20110811.root","RECREATE");
 TTree *tout=new TTree("tree_passedEvents","Converted from FP");
 tout->Branch("CMS_hzz2l2q_mZZ",&mZZOUT,"CMS_hzz2l2q_ZZ/D");
 tout->Branch("mZjj",&mZjjOUT,"mZjj/D");
 tout->Branch("leptType",&leptTypeOUT,"leptType/D");
 tout->Branch("nBTags",&nBTagsOUT,"nBTags/D");
 tout->Branch("eventWeight",&eventWeightOUT,"eventWeight/D");

 cout<<"start"<<endl;
 for(int i=0;i<t->GetEntries();i++){
   t->GetEntry(i);
   leptTypeOUT=double(leptTypeIN);
   nBTagsOUT=double(nBTagsIN);
   mZZOUT=double(mZZIN);
   mZjjOUT=double(mZjjIN);
   eventWeightOUT=double(eventWeightIN);
   tout->Fill();
 }
 cout<<"finished"<<endl;
 tout->Write();
 delete tout;
 delete fout;

}
