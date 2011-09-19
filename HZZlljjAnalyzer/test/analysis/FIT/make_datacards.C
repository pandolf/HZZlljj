#include <Riostream.h>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

#include "TGraph.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"

#include "make_roofitfiles.C"
#include "HiggsCSandWidth.cc"

const int nbtag=3; //btag0, btag1 and btag2
const int nprod=2; //VBF and gg
const int nchan=2;//2e2j and 2m2j
const float lumiee=1.556; //fb^-1, NEW!
const float lumimm=1.615; //fb^-1, NEW!
const int nmass=15;
const bool isSM4=false; // if true add .15 to CSgg and CSvbf errors 
const float mass[nmass]={190.0,200.0,210.0,230.0,250.0, 
                         300.0,350.0,400.0,425.0,450.0,
                         500.0,525.0,550.0,575.0,600.0};


float exp_sig_yields[nbtag][nchan][nmass];
//float exp_sig_yields_ee[nbtag][nmass];
//float exp_sig_yields_mm[nbtag][nmass];

float gen_sig_yields[nmass];
//={//no mass cut here
//  264741, 284748, 291226, 109981,
//  296279, 109076, 109924, 291650, 109906,
//  102857, 286612, 292613, 299930, 109977
//};

const int maxNbkg[nbtag]={1,1,1};//for every btag, there are MAX==2 different types of bkg (only one for btag0 and 1, ttbar +all the rest for btag2)
float obs_bkg[nbtag][nchan]={//from real data
  307, 359,//btag0, first 2e2j, then 2m2j
  351, 396,
  30,  25
 };

float exp_bkg[nbtag][nchan];
//={//from sideband procedure 
//  286.4, 345.7,//btag0, first 2e2j, then 2m2j
//  334.7, 376.4,//btag1
//  20.3,  24.3// 15.3,  29.7  //btag2
//};
 
//const float NLOxsect[nmass]={ 144.1,106.3, 80.2, 76.9, 55.3, 35.3, 22.4, 8.4};//in fb
//const float meas_lumi[nmass]={ 757.2,1034.6, 1360.4, 1418.4, 1987.8, 3113.5, 4591.8,109977/ 8.4};//in fb^-1
vector<float> btag0_ee_eff_pars,  btag1_ee_eff_pars, btag2_ee_eff_pars; //we store the parameters of the dependence of the eff on the Mzz
vector<float> btag0_mm_eff_pars,  btag1_mm_eff_pars, btag2_mm_eff_pars; //we store the parameters of the dependence of the eff on the Mzz
const float maxDiffSyst=0.01;
HiggsCSandWidth *myCSW;

//vector<float>aux;


string make_ratestring(float mymass, vector<float> myxsect, int ibtag,int chan, float expbkg=-99.9);
string make_theorunc(float mymass,float myerrp, float myerrm,bool is_ggH=true, int ibtag=0);  
string make_theorunc(float mymass,float myerrCSp, float myerrCSm, float myerrVBFp, float myerrVBFm);
string make_pdftheorunc(float mymass,float myerrp, float myerrm,bool is_ggH=true, int ibtag=0);  
string make_theorgammaunc(float mymass);
string make_btagunc(float mymass,int ibtag, int ich); 
string make_JESunc(float mymass);
string make_bakcgrNormErrLine(double expyields, int ibtag, int ich); 
string make_obsstring(int btag);
string make_obsstring(int btag, int obs);
vector<float> eff_fit(int mybtag,int chan);
vector<double> calculate_CBpars(int btag,double mH);
int get_mass_index(float mymass);
void extract_exp_sig_yields(float mymass,float mywidth);//,vector<float> ee_y , vector<float> mm_y);
void get_gen_yields();
int make_datacards();
int get_signal_gen(float mymass);
int make_datacards(){

  // gROOT->ProcessLine(".L make_roofitfiles.C++");

  myCSW = new HiggsCSandWidth();

  // return 1;
 
  string nameXsecFile;
  if(isSM4)
    nameXsecFile="./xsect_higgs_173points_4generation.txt";
  else
    nameXsecFile="./xsect_higgs_173points_new.txt";
  ifstream xsect_file(nameXsecFile.c_str(),ios::in);

  if (! xsect_file.is_open()){ cout<<"Failed to open file with xsections"<<endl;}
  float mH, CSgg, CSgg_p, CSgg_m, CSpdfgg_p,CSpdfgg_m,CSvbf, CSvbf_p, CSvbf_m,CSpdfvbf_p,CSpdfvbf_m, 
        Gamma, BRHZZ, BRZZ2l2q;

  string outDir="datacards/";

  //calculate expected yields
  get_gen_yields();
  while(xsect_file.good()){
    xsect_file >> mH>> CSgg>> CSgg_p >> CSgg_m >>  CSpdfgg_p>>CSpdfgg_m>>CSvbf >> CSvbf_p >> CSvbf_m >>  CSpdfvbf_p>>CSpdfvbf_m>>Gamma >> BRHZZ >> BRZZ2l2q;

    //    cout<<mH<<"  "<<CSgg<<endl;
    //if(mH<250.0) continue;
    if(mH==420)mH=425.;
    if(mH==460)mH=450.;
    if(mH==520)mH=525.;
    if(mH==560)mH=550.;
    if(mH==580)mH=575.;
    extract_exp_sig_yields(mH,Gamma);//exp_sig_yields_ee,exp_sig_yields_mm);
    extract_exp_sig_yields(mH,Gamma);//exp_sig_yields_ee,exp_sig_yields_mm);
 }

  //parametrize efficiencies
  btag0_ee_eff_pars=eff_fit(0,0);
  //  cout<<"In main: "<<btag0_eff_pars.at(0)<<"    "<<btag0_eff_pars.at(1)<<endl;
  btag1_ee_eff_pars=eff_fit(1,0);
  btag2_ee_eff_pars=eff_fit(2,0);
  btag0_mm_eff_pars=eff_fit(0,1);
  btag1_mm_eff_pars=eff_fit(1,1);
  btag2_mm_eff_pars=eff_fit(2,1);
 
  xsect_file.clear();
  xsect_file.seekg (0);

  while(xsect_file.good()){
    xsect_file >> mH>> CSgg>> CSgg_p >> CSgg_m >>  CSpdfgg_p>>CSpdfgg_m>>CSvbf >> CSvbf_p >> CSvbf_m >>  CSpdfvbf_p>>CSpdfvbf_m>>Gamma >> BRHZZ >> BRZZ2l2q;
    //    cout<<mH<<"  "<<CSgg<<endl;

    // ========   add additional uncertainty to cs for SM4
    if(isSM4){
      CSgg_p+=.1;
      CSgg_m-=.1;
      CSvbf_p+=.1;
      CSvbf_m-=.1;
    }// modified by A Whitbeck

    if(mH<200.0) continue;
    // if(mH!=300.0) continue;
    //create directory for this mass point
    std::ostringstream ossDir;
    ossDir<<mH;
    if(mkdir((outDir+ossDir.str()).c_str(),0777)==-1){
    //if(mkdir((outDir+ossDir.str()).c_str())==-1){
      //cout<<"Failed to create directory "<<ossDir.str()<<endl;
      //cout<<"Failed to create directory "<<outDir+ossDir.str().c_str()<<endl;
      //break;
    }

    float CSggxs =1000.0*CSgg  *BRHZZ *BRZZ2l2q;
    float CSvbfxs=1000.0*CSvbf*BRHZZ  *BRZZ2l2q;
    //   cout<<"THEOR xsect: "<<CSgg<<" , Vbf= "<<CSvbf<<"  BRHZZ="<<BRHZZ<<"  BRZZ2L2q="<<BRZZ2l2q<<endl;
    vector<float> myxsect;
    myxsect.push_back(CSggxs);
    myxsect.push_back(CSvbfxs);

    float ggRelErrp=CSgg_p;
    float ggRelErrm=CSgg_m;//the error quoted in the txt file has already a minus
    float vbfRelErrp=CSvbf_p;
    float vbfRelErrm=CSvbf_m;
    float ggPDFErrp=CSpdfgg_p;
    float ggPDFErrm=CSpdfgg_m;//the error quoted in the txt file has already a minus
    float vbfPDFErrp=CSpdfvbf_p;
    float vbfPDFErrm=CSpdfvbf_m;
  //cout<<"MH = "<<mH<<"    "<<test1.c_str()<<endl;  

    string str_ch[2]={"ee","mm"};
    string str_btag[3]={"0b","1b","2b"};
      //loop over btag categories
    for(int ibtag=0;ibtag<nbtag;ibtag++){
   
      string theorErrLine_ggH=make_theorunc(mH, ggRelErrp, ggRelErrm, true,ibtag);  
      string theorErrLine_VBF=make_theorunc(mH, vbfRelErrp, vbfRelErrm, false,ibtag);  
      string pdfErrLine_ggH=make_pdftheorunc(mH, ggPDFErrp, ggPDFErrm, true,ibtag); 
      string pdfErrLine_VBF=make_pdftheorunc(mH, vbfPDFErrp, vbfPDFErrm, false,ibtag);  
      string gammaErrLine=make_theorgammaunc(mH);

      string jesLine=make_JESunc(mH);
      //loop over decay channels
      for(int ich=0;ich<nchan;ich++){
 
	string beffLine=make_btagunc(mH,ibtag,ich);
  	string str_id=str_ch[ich]+str_btag[ibtag];

      ifstream tpl_file(("./hzz2l2q_"+ str_id+".tpl").c_str(),ios::in);

	std::ostringstream mass_str;
	mass_str<<mH;

	//ok, now produce the root file wit hthe workspace
	vector<double> cbpars =calculate_CBpars(ibtag,mH);
	double expyields=0.0;// exp_bkg[ibtag][ich];
	double obsyields=0.0;
	make_roofitfiles(ibtag,ich,mH,Gamma, obsyields, expyields,cbpars);
	exp_bkg[ibtag][ich]=expyields;
	cout<<"\n\n~~~~~~~~~~Btag "<< ibtag<<"   Expected MAKEDATACARDS "<<expyields<<endl;
	string expRateLine=make_ratestring(mH,myxsect,ibtag,ich,expyields);
	string obsLine=make_obsstring(ibtag, int(obsyields));
	string bakcgrNormErrLine = make_bakcgrNormErrLine(expyields,ibtag,ich);

	string datacard_name=outDir+ossDir.str()+"/hzz2l2q_"+str_id+"."+mass_str.str() +".txt";
	ofstream datacard_new(datacard_name.c_str(),ios::out);
	string tpl_line;
	string str_ratetag("<dummy1>");
	string str_obstag("<dummyobs>");
	string str_befftag("<dummybeff>");
	string str_gammatag("<dummygammaBW>");
	string str_jestag("<dummyJES>");
	string str_theorggHtag("<dummyggH>");
	string str_theorVBFtag("<dummyVBF>");
	string str_pdftheorggHtag("<dummypdfggH>");
	string str_pdftheorVBFtag("<dummypdfqqH>");
	string str_bakcgrNorm("<dummybnorm>");
	string str_bckgShape_0b_1("CMS_hzz2l2q_bkg0bp1");
	string str_bckgShape_0b_2("CMS_hzz2l2q_bkg0bp2");
	string str_bckgShape_0b_3("CMS_hzz2l2q_bkg0bp3");
	string str_bckgShape_0b_4("CMS_hzz2l2q_bkg0bp4");
	string str_bckgShape_1b_1("CMS_hzz2l2q_bkg1bp1");
	string str_bckgShape_1b_2("CMS_hzz2l2q_bkg1bp2");
	string str_bckgShape_1b_3("CMS_hzz2l2q_bkg1bp3");
	string str_bckgShape_1b_4("CMS_hzz2l2q_bkg1bp4");
	string str_bckgShape_2b_1("CMS_hzz2l2q_bkg2bp1");
	string str_bckgShape_2b_2("CMS_hzz2l2q_bkg2bp2");
	string str_bckgShape_2b_3("CMS_hzz2l2q_bkg2bp3");
	string str_bckgShape_2b_4("CMS_hzz2l2q_bkg2bp4");
	bool found=false;
	while(tpl_file.good()){
	  getline (tpl_file,tpl_line);
	  size_t posrate_found=tpl_line.find(str_ratetag);
	  size_t posobs_found=tpl_line.find(str_obstag);
	  size_t posbeff_found=tpl_line.find(str_befftag);
	  size_t posjes_found=tpl_line.find(str_jestag);
	  size_t postheorgamma_found=tpl_line.find(str_gammatag);
	  size_t postheorggH_found=tpl_line.find(str_theorggHtag);
	  size_t postheorVBF_found=tpl_line.find(str_theorVBFtag);
	  size_t pospdfggH_found=tpl_line.find(str_pdftheorggHtag);
	  size_t pospdfVBF_found=tpl_line.find(str_pdftheorVBFtag);
	  size_t posbnorm_found=tpl_line.find(str_bakcgrNorm);
	  size_t pos_bckgShape_0b_1_found=tpl_line.find(str_bckgShape_0b_1);
	  size_t pos_bckgShape_0b_2_found=tpl_line.find(str_bckgShape_0b_2);
	  size_t pos_bckgShape_0b_3_found=tpl_line.find(str_bckgShape_0b_3);
	  size_t pos_bckgShape_0b_4_found=tpl_line.find(str_bckgShape_0b_4);
	  size_t pos_bckgShape_1b_1_found=tpl_line.find(str_bckgShape_1b_1);
	  size_t pos_bckgShape_1b_2_found=tpl_line.find(str_bckgShape_1b_2);
	  size_t pos_bckgShape_1b_3_found=tpl_line.find(str_bckgShape_1b_3);
	  size_t pos_bckgShape_1b_4_found=tpl_line.find(str_bckgShape_1b_4);
	  size_t pos_bckgShape_2b_1_found=tpl_line.find(str_bckgShape_2b_1);
	  size_t pos_bckgShape_2b_2_found=tpl_line.find(str_bckgShape_2b_2);
	  size_t pos_bckgShape_2b_3_found=tpl_line.find(str_bckgShape_2b_3);
	  size_t pos_bckgShape_2b_4_found=tpl_line.find(str_bckgShape_2b_4);
	  if(posrate_found!=string::npos){
	    //	 cout<<tpl_line.c_str()<<endl;
	    found=true;
	    datacard_new<<expRateLine<<endl;
	  }
	  else if(posobs_found!=string::npos){
	    datacard_new<<obsLine<<endl;
	  }
	  else if(posbeff_found!=string::npos){
	    datacard_new<<beffLine<<endl;
	  }
	  else if(posjes_found!=string::npos){
	    datacard_new<<jesLine<<endl;
	  }
	  else if(postheorggH_found!=string::npos){
	    datacard_new<<theorErrLine_ggH<<endl;
	  }
	  else if(postheorVBF_found!=string::npos){
	    datacard_new<<theorErrLine_VBF<<endl;
	  }
	  else if(pospdfggH_found!=string::npos){
	    datacard_new<<pdfErrLine_ggH<<endl;
	  }
	  else if(pospdfVBF_found!=string::npos){
	    datacard_new<<pdfErrLine_VBF<<endl;
	  }
	   else if(posbnorm_found!=string::npos){
	    datacard_new<<bakcgrNormErrLine<<endl;
	   }
	  /*else if(pos_bckgShape_0b_1_found!=string::npos  || pos_bckgShape_0b_2_found!=string::npos ||
		  pos_bckgShape_0b_3_found!=string::npos  || pos_bckgShape_0b_4_found!=string::npos ||
		  pos_bckgShape_1b_1_found!=string::npos  || pos_bckgShape_1b_2_found!=string::npos ||
		  pos_bckgShape_1b_3_found!=string::npos  || pos_bckgShape_1b_4_found!=string::npos ||
		  pos_bckgShape_2b_1_found!=string::npos  || pos_bckgShape_2b_2_found!=string::npos ||
		  pos_bckgShape_2b_3_found!=string::npos  || pos_bckgShape_2b_4_found!=string::npos ){
	    datacard_new<<"# no shape errors for now"<<endl;
	    }*/
	  //else if(postheorgamma_found!=string::npos){
	  //   datacard_new<<gammaErrLine<<endl;
	  //  }
	  else datacard_new<<tpl_line<<endl;
	}//ed loop over lines of template file
	if(!found)cout<<"it did not find the tag string "<<str_ratetag.c_str()<<endl;
	


      }//end loop on channels
    }//end loop on btag
  }//end loop over xsect file (-> loop over mass points)
  

  
  return 0;
}//end main
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////

string make_ratestring(float mymass,vector<float> myxsect, int ibtag, int chan,float expbkg){

  //vector of signal eff: btag0_2e2j_VBF,  btag0_2e2j_gg,  btag0_2m2j_VBF,  btag0_2m2j_gg,  btag1_2e2j_VBF, ...  btag2_2m2j_gg = 12 entries
  vector<float> sig_eff;
  string rate_str;
  rate_str="rate     ";
  //signal efficiency parametrization
  float lumi=0.0;
  float eff=0.0;
  float mya=0.0, myb=0.0, myc=0.0,myd=0.0;
  if(chan==0){
    if(ibtag==0){
      mya=btag0_ee_eff_pars.at(0);
      myb=btag0_ee_eff_pars.at(1);
      myc=btag0_ee_eff_pars.at(2);
      myd=btag0_ee_eff_pars.at(3);
    }
    else if(ibtag==1){
      mya=btag1_ee_eff_pars.at(0);
      myb=btag1_ee_eff_pars.at(1);
      myc=btag1_ee_eff_pars.at(2);
      myd=btag1_ee_eff_pars.at(3);
    }
    else if(ibtag==2){
      mya=btag2_ee_eff_pars.at(0);
      myb=btag2_ee_eff_pars.at(1);
      myc=btag2_ee_eff_pars.at(2);
      myd=btag2_ee_eff_pars.at(3);
    }
    else{
      cout<<"\n\nERROR !! Btag index out of range!"<<endl<<endl<<endl;
      return "dummy";
    }
    lumi=lumiee;
  }
  else{//chan==1 -> muons
    if(ibtag==0){
      mya=btag0_mm_eff_pars.at(0);
      myb=btag0_mm_eff_pars.at(1);
      myc=btag0_mm_eff_pars.at(2);
      myd=btag0_mm_eff_pars.at(3);
    }
    else if(ibtag==1){
      mya=btag1_mm_eff_pars.at(0);
      myb=btag1_mm_eff_pars.at(1);
      myc=btag1_mm_eff_pars.at(2);
      myd=btag1_mm_eff_pars.at(3);
    }
    else if(ibtag==2){
      mya=btag2_mm_eff_pars.at(0);
      myb=btag2_mm_eff_pars.at(1);
      myc=btag2_mm_eff_pars.at(2);
      myd=btag2_mm_eff_pars.at(3);
    }
    else{
      cout<<"\n\nERROR !! Btag index out of range!"<<endl<<endl<<endl;
      return "dummy";
    }
    lumi=lumimm;
  }
  //  eff=mya + myb*mymass;
  //cout<<"BTAG = "<<ibtag<<" CHAN = "<<chan<<" using these params: "<<mya<<" "<<myb<<"  "<<myc<<"  "<<myd<<"  Mass: "<<mymass<<endl;
  eff=mya + myb*mymass + myc *mymass*mymass + myd *mymass*mymass*mymass;
  //  eff = mya+myb*log(mymass);
  //read from text file the right Higgs xsect
  //loop over prod mechanism
  for(int iprod=0;iprod<nprod;iprod++){
    //    cout<<"\n\nBTAG "<<ibtag<<"  CHAN "<<chan<<"  PROD "<<iprod<<"   :  Eff= "<<eff<<endl;
    float xsect=myxsect.at(iprod);
   

    // cout<<"\n\n~~~M = "<<mymass<<"BTAG = "<<ibtag<<" CHAN = "<<chan<<"  eff="<<eff<<"  xsect="<<xsect<<"  lumi="<<lumi<<"  BrZtoL="<<brZtoL<<flush;
    float nexpected=eff*xsect*lumi*0.5;//divide by two because the xsect is the sum of 2e2j and 2m2j
    //cout<<"  ---> EXP: "<<nexpected<<endl;
    sig_eff.push_back(nexpected);
    std::ostringstream ossm;
    ossm<<nexpected;
    string nexp_str=ossm.str() ;
    rate_str+=("    "+nexp_str);
    if(iprod==nprod-1){//finished to fill signal yields, add bakgd
      
      // for(int ibkg=0;ibkg<maxNbkg[ibtag];ibkg++){//not necessary any more, all bkg lumped into one also for 2btag
	std::ostringstream ossbkg;
	if(expbkg<0.0)ossbkg<<exp_bkg[ibtag][chan];
	else ossbkg<<expbkg;
	string bkg_str=ossbkg.str() ;
	rate_str+=("    "+bkg_str);
	//	  }//end loop on diferent types of bkgd
    }//end if iprod==nprod-1
    
    //close three loops
      }//end loop on prod mechanisms
    
  //background yields are fixed (because we fit the shape over the whole mass range)
  
  return rate_str;
}//END string make_ratestring(float mass)



string make_theorunc(float mymass, float myerrp, float myerrm, bool is_ggH, int ibtag){
  //ERRORS SHOULD BE RELATIVE AND NETTO, I.E.: for a 6% syst, pass 0.06 AND NOT 1.06 
  string theor_str;
  if(is_ggH)  theor_str="QCDscale_ggH    lnN       ";
  else theor_str="QCDscale_qqH    lnN       "; //else is VBF
  //float therrp=0., therrm=0.0;
  string nthp="-999";
  string nthm="+999";
  std::ostringstream ossp;
  myerrp+=1.0;
  myerrm+=1.0;
  ossp<<(myerrp);
  nthp=ossp.str() ;
  std::ostringstream ossm;
  ossm<<(myerrm);//negative errors are passed already with the minus sign
  nthm=ossm.str() ;
  string therr_str=nthm+"/"+nthp;
  bool asymm=true;
  float erravg=(fabs(myerrp-1.0)+fabs(myerrm-1.0))/2.0 + 1.0;
  if(fabs(erravg-myerrp)/erravg <maxDiffSyst ){
    asymm=false;//symmetric errors
    std::ostringstream ossavg;
    ossavg<<erravg;
    therr_str=ossavg.str();
  }

  if(is_ggH){
    theor_str+=("    "+therr_str);
    //for the VBF part, no error
    theor_str+="    1.0";
    //for the background, no err
    theor_str+="    1.0"; 
    //  if(ibtag==2){//two bakg sources for the btag2 category
    //    theor_str+="    1.0";
    //  }
  }
  else{
    //for the ggH part, no error
    theor_str+="    1.0";
    theor_str+=("    "+therr_str);
    //for the background, no err
    theor_str+="    1.0"; 
    //  if(ibtag==2){//two bakg sources for the btag2 category
      //    theor_str+="    1.0";
    //  }
  }
  //  if(int(mymass)%50==0) cout<<"*** MASS "<<mymass<<" "<<theor_str.c_str()<<endl;
  return theor_str;
}//end make_theorunc
  

string make_theorunc(float mymass,float myerrCSp, float myerrCSm, float myerrVBFp, float myerrVBFm){
 //ERRORS SHOULD BE RELATIVE AND NETTO, I.E.: for a 6% syst, pass 0.06 AND NOT 1.06 
  string theor_str;
  theor_str="QCDscale_ggH    lnN       ";
 //loop over btag categories
  for(int ibtag=0;ibtag<nbtag;ibtag++){
    //loop over decay channels
    for(int ich=0;ich<nchan;ich++){
      //float therrp=0., therrm=0.0;
     //loop over prod mechanism
      for(int iprod=0;iprod<nprod;iprod++){// same err for all chans and btags
	string nthp="-999";
	string nthm="+999";
	if(iprod==0){//gg	  
	  std::ostringstream ossp;
	  myerrCSp+=1.0;
	  myerrCSm+=1.0;
	  ossp<<(myerrCSp);
	  nthp=ossp.str() ;
	  std::ostringstream ossm;
	  ossm<<(myerrCSm);//negative errors are passed already with the minus sign
	  nthm=ossm.str() ;
	}
	else if(iprod==1){//vbf	  
	  std::ostringstream ossp;
	  myerrVBFp+=1.0;
	  myerrVBFm+=1.0;
	  ossp<<(myerrVBFp);
	  nthp=ossp.str() ;
	  std::ostringstream ossm;
	  ossm<<(myerrVBFm);
	  nthm=ossm.str() ;
	}
	else{
	  cout<<"IPROD out of range"<<endl;
	}
	theor_str+=("    "+nthm+"/"+nthp);
	//for the background, no err 
	if(iprod==nprod-1){
	  theor_str+="    1.0";
	  // if(ibtag==2){//two bakg sources for the btag2 category
	  //    theor_str+="    1.0";
	  //  }
	}
      }//end  for(int iprod
    }//end  for(int ich=0;
  }//end  for(int ibtag

  return theor_str;
}//end make_theorunc

string make_pdftheorunc(float mymass, float myerrp, float myerrm, bool is_ggH, int ibtag){
 //ERRORS SHOULD BE RELATIVE AND NETTO, I.E.: for a 6% syst, pass 0.06 AND NOT 1.06 
  string theor_str;
  if(is_ggH)  theor_str="pdf_gg   lnN       ";
  else theor_str="pdf_qqbar    lnN       "; //else is VBF
  //float therrp=0., therrm=0.0;
  string nthp="-999";
  string nthm="+999";
  float signp = (myerrp>0 ? 1.0 : -1.0 );
  float signm = (myerrm>0 ? 1.0 : -1.0 );

  float add_vbf_err=0.0;
  /*if(!is_ggH){//for VBF, add errors due to imperfect knowledge of eff, see Table 23 of AN-2011/100
   
    if(ibtag==0){
      if(mymass<=500.0)add_vbf_err=0.01;
      else add_vbf_err=0.02;
    }
    else if(ibtag==1){
      if(mymass<=500.0)add_vbf_err=0.01;
      else add_vbf_err=0.04;
    }
    else if(ibtag==2){
      if(mymass<=500.0)add_vbf_err=0.01;
      else add_vbf_err=0.02;
    }
    else{
      cout<<"wrong btag passed as argument to make_pdftheorunc"<<endl;
    }
    }*/ //has nothing to do with PDF!!!

  std::ostringstream ossp;
  float errplus=signp*sqrt(myerrp*myerrp+0.015*0.015+add_vbf_err*add_vbf_err);//0.015 is the syst on signal acceptance due to pdf uncertainties
  float errminus=signm*sqrt(myerrm*myerrm+0.04*0.04+add_vbf_err*add_vbf_err);//0.04 is the syst on signal acceptance due to pdf uncertainties
  errplus=errplus+1.0;
  errminus=errminus+1.0;
  ossp<<errplus;
  nthp=ossp.str() ;
  std::ostringstream ossm;
  ossm<<errminus;//negative errors are passed already with the minus sign
  nthm=ossm.str() ;
  string therr_str=nthm+"/"+nthp;
  bool asymm=true;
  float erravg=(fabs(errplus-1.0)+fabs(errminus-1.0))/2.0 +1.0;
   if(fabs(erravg-errplus)/erravg <maxDiffSyst ){
    asymm=false;//symmetric errors
    std::ostringstream ossavg;
    ossavg<<erravg;
    therr_str=ossavg.str();
     }

  if(is_ggH){
    theor_str+=("    "+therr_str);
    //for the VBF part, no error
    theor_str+="    1.0";
    //for the background, no err
    theor_str+="    1.0"; 
    // if(ibtag==2){
    //  theor_str+="    1.0";
    // }
  }
  else{
    //for the ggH part, no error
    theor_str+="    1.0";
    theor_str+=("    "+therr_str);
    //for the background, no err
    theor_str+="    1.0"; 
    // if(ibtag==2){//two bakg sources for the btag2 category
    //   theor_str+="    1.0";
    // }
  }

  //  cout<<"+++ MASS "<<mymass<<" "<<theor_str.c_str()<<endl;
  return theor_str;
}//end make_pdftheorunc

string make_theorgammaunc(float mymass){
  string gamma_str="theory_gamma  lnN  ";

  float mHinTeV=mymass/1000.0;
  float unc=150*(mHinTeV*mHinTeV*mHinTeV);//percentual error on xsect
  std::ostringstream ossp;
  ossp<<(1.0+unc/100.0);
  string tmp_str=ossp.str()+"   "+ossp.str()+"   1.0" ;
  gamma_str+=tmp_str;
  return gamma_str;
}


string make_obsstring(int btag){//superseded !!!!!
  string obs_str="observation   ";

  std::ostringstream ossp;
  ossp<<int(exp_bkg[btag][0]/2.0);
  obs_str+=ossp.str() ;

  return obs_str;
}

string make_obsstring(int btag, int obs){
  string obs_str="observation   ";

  std::ostringstream ossp;
  ossp<<obs;
  obs_str+=ossp.str() ;

  return obs_str;
}

string make_btagunc(float mymass,int ibtag, int ich){
  string b_str="CMS_eff_b	lnN";
  float p0=0.0, p1=0.0;
  float m0=0.0, m1=0.0;
  float errp=999.0, errm=-9999.0;
  if(ich==0){ //electrons
    if(ibtag==0){
      p0=0.983256647923;
      p1=-0.0000883532570978;
      m0=1.02907356239;
      m1=0.0000713061639147;
    }
    else if(ibtag==1){
      p0=1.04446110956;
      p1=-0.0000195508160829;
      m0=0.940063743731;
      m1=0.0000737044467898;
    }
    else if(ibtag==2){
      p0=1.13365470372;
      p1=0.00000584572717646;
      m0=0.82161771535;
      m1=-0.0000161054152592;
    }
    else{
      cout<<"Wrong # btags to make_btagunc"<<endl;
    }
  }
  else if(ich==1){ //muons
    if(ibtag==0){
      p0=0.984636818312;
      p1=-0.0000898705296203;
      m0=1.02836905579;
      m1= 0.0000726807344479;
    }
    else if(ibtag==1){
      p0=1.04385580002;
      p1=0.0000206096278947;
      m0=0.942713582987;
      m1=0.0000719882385098;
    }
    else if(ibtag==2){
      p0=1.1333366687;
      p1=0.00000462542786413;
      m0=0.813316607701;
      m1=-0.00000205840248842;
    }
    else{
      cout<<"Wrong # btags to make_btagunc"<<endl;
    }
  }
  else{
      cout<<"Wrong # channel to make_btagunc"<<endl;
    }
  cout<<"p0 "<<p0<<endl;
  cout<<"p1 "<<p1<<endl;
  cout<<"m0 "<<m0<<endl;
  cout<<"m1 "<<m1<<endl;
  cout<<"mymass "<<mymass<<endl;
  errp=p1*mymass+p0;
  errm=m1*mymass+m0;
  cout<<"errp "<<errp<<endl;
  cout<<"errm "<<errm<<endl;
  float erravg=(fabs(errp-1.0)+fabs(errm-1.0))/2.0 +1.0;
  cout<<"erravg "<<erravg<<endl;
  bool asymm=true;
  if(fabs(erravg-errp)/erravg <maxDiffSyst ){
    asymm=false;//symmetric errors
    // cout<<"@@@@@@@ asymm errors for BTAG ("<< maxDiffSyst<<"): "<<errm<<"  "<<errp<<"  "<<erravg<<endl;
  }

  std::ostringstream ossp;
  ossp<<errp;
  std::ostringstream ossm;
  ossm<<errm;
  std::ostringstream ossavg;
  ossavg<<erravg;
  cout<<"errp "<<errp<<endl;
  cout<<"errm "<<errm<<endl;
  cout<<"erravg "<<erravg<<endl;
  string tmp_str;//=ossm.str()+"/"+ossp.str();
  if(asymm)tmp_str=ossm.str()+"/"+ossp.str();
  else tmp_str=ossavg.str();
  b_str+=("  "+tmp_str+"      "+tmp_str+"      1.0");
  return b_str;
}

string make_JESunc(float mymass){

  string jes_str="CMS_scale_j	lnN";
  float p0= 8.3  , p1=-0.0215 ;
  float m0=-8.6, m1=0.02 ;
  float errp=999.0, errm=-9999.0;

  //in percent
  errp=1.0+0.01*(p0+p1*mymass);
  errm=1.0+0.01*(m0+m1*mymass);
  //  cout<<"========= ERR JES:  "<<errp<<"   "<<errm<<endl;
  float erravg=(errp-1.0+fabs(errm-1.0))/2.0 +1.0;
  bool asymm=true;
  if(fabs(erravg-errp)/erravg <maxDiffSyst ){
    asymm=false;//symmetric errors
  }

  std::ostringstream ossp;
  ossp<<errp;
  std::ostringstream ossm;
  ossm<<errm;
  std::ostringstream ossavg;
  ossavg<<erravg;
  string tmp_str;//=ossm.str()+"/"+ossp.str();
  if(asymm)tmp_str=ossm.str()+"/"+ossp.str();
  else tmp_str=ossavg.str();
  jes_str+=("  "+tmp_str+"      "+tmp_str+"      1.0");
  return jes_str;
}

string make_bakcgrNormErrLine(double expyields, int btag, int ich){
  string str_btag,str_chann,str_neventsSB, str_alpha;
  double alpha;
  if(ich==0)
    str_chann="ee";
  else if(ich==1)
    str_chann="mm";
  if(btag==0){
    str_btag="0b";
    if(ich==0){ //electrons
      str_neventsSB="529";
      alpha=expyields/529.;
      std::ostringstream ossa;
      ossa<<alpha;
      str_alpha=ossa.str();
    }
    else if(ich==1){ //muons
      str_neventsSB="621";
      alpha=expyields/621.;
      std::ostringstream ossa;
      ossa<<alpha;
      str_alpha=ossa.str();
    }
  }
  else if(btag==1){
    str_btag="1b";
    if(ich==0) {//electrons
      str_neventsSB="534";
      alpha=expyields/534.;
      std::ostringstream ossa;
      ossa<<alpha;
      str_alpha=ossa.str();
    }
    else if(ich==1){ //muons
      str_neventsSB="615";
      alpha=expyields/615.;
      std::ostringstream ossa;
      ossa<<alpha;
      str_alpha=ossa.str();
    }
  }
  else if(btag==2){
    str_btag="2b";
    if(ich==0) {//electrons
      str_neventsSB="37";
      alpha=expyields/37.;
      std::ostringstream ossa;
      ossa<<alpha;
      str_alpha=ossa.str();
    }
    else if(ich==1) {//muons
      str_neventsSB="44"; 
      alpha=expyields/44.;
      std::ostringstream ossa;
      ossa<<alpha;
      str_alpha=ossa.str();
    }  
  }
  string bNorm_str="CMS_hzz2l2q_bkg"+str_btag+str_chann+"p0    gmN   "+str_neventsSB+" ----  -----  "+str_alpha;
  return bNorm_str;
}

vector<float> eff_fit(int mybtag,int chan){
  cout<<"\n~~~~~~~ Eff fit for btag="<<mybtag<<" chan = "<<chan<<endl;
  std::ostringstream ossm;
  ossm<<mybtag;
  string btag_str="btag"+ossm.str();
  string chan_str= (chan==0)?"ee" : "mm";
  float myeff[nmass];
  for(int im=0;im<nmass;im++){
    //  myeff[im]=exp_eff[mybtag][im];

    float brZtoL=0.0;
    if(get_signal_gen(mass[im])==0){//for powheg samples, Z->ll with l=e,mu,tau; in jhu, Z->ll with l=e,mu
      brZtoL=1.0/3.0;
    }
    else brZtoL=0.5;
   
    myeff[im]=exp_sig_yields[mybtag][chan][im]/(gen_sig_yields[im]*brZtoL);  //(lumi*NLOxsect[im]);
    //  cout<<"Calc EFF: im="<<im<<"  chan="<<chan<<" "<<mybtag<<"b  EXP="<<exp_sig_yields[mybtag][chan][im]<<"  GEN="<<gen_sig_yields[im]<<"  brZtoL="<<brZtoL<<" ---> EFF: "<<myeff[im]<<endl;
    //if(chan==0)  myeff[im]=exp_sig_yields[mybtag][chan][im]/(gen_sig_yields[im]/2.0);  //(lumi*NLOxsect[im]);
    //else  myeff[im]=exp_sig_yields_mm[mybtag][im]/(gen_sig_yields[im]/2.0);//(lumi*NLOxsect[im]);
    // cout<<"Eff for mass "<<mass[im]<<"  BTag "<<mybtag <<"  "<<myeff[im]<<endl;
  }

  TGraph *gr_eff=new TGraph(nmass,mass,myeff);
  gr_eff->SetName(("effgr_"+btag_str+chan_str).c_str());
  gr_eff->SetTitle(("Efficiency vs Mass ("+chan_str+" , "+btag_str+")").c_str());
  gr_eff->SetMarkerStyle(21);
  gr_eff->SetMinimum(0.0);

  // TF1 *f1 = new TF1("fit_pol1","pol1",1,2);
  // TF1 *f1 = new TF1("fit_sqrt","[0]+[1]*sqrt(x)",1,2);
  // TF1 *f1 = new TF1("fit_log","[0]+[1]*log(x)",150.0,800.0);
  string f1name="fit_poly3_"+btag_str+chan_str;
  TF1 *f1 = new TF1(f1name.c_str(),"pol3",150.0,700.0);
  gr_eff->Fit(f1name.c_str(),"Q");
  vector<float> fit_res;
  fit_res.push_back(f1->GetParameter(0) );
  fit_res.push_back(f1->GetParameter(1) );
  fit_res.push_back(f1->GetParameter(2) );
  fit_res.push_back(f1->GetParameter(3) );
  cout<<"\n\nFit results for BTAG "<<mybtag<<"  (eff= a + b*Mzz + c*Mzz^2 + d*Mzz^3):"<<endl;
  cout<<" a = "<<fit_res.at(0)<<" +/- "<<f1->GetParError(0) <<endl;
  cout<<" b = "<<fit_res.at(1)<<" +/- "<<f1->GetParError(1) <<endl;
  cout<<" c = "<<fit_res.at(2)<<" +/- "<<f1->GetParError(2) <<endl;
  cout<<" d = "<<fit_res.at(3)<<" +/- "<<f1->GetParError(3) <<endl<<endl;

  string canname="c_fitfunc_"+btag_str+chan_str;
  TCanvas *c1=new TCanvas(canname.c_str(),"CANFIT", 900,900);
  c1->SetFillColor(kWhite);
  c1->SetBorderMode(0);
  c1->cd();
  gPad->SetFillColor(kWhite);
  gr_eff->GetXaxis()->SetTitle("M_{H}");
  gr_eff->GetYaxis()->SetTitle("#varepsilon");
  gr_eff->GetYaxis()->SetTitleOffset(1.0);
  gr_eff->Draw("AP");
  c1->SaveAs( ("eff_param_newHLT_"+btag_str +chan_str+".root").c_str()  );

  delete c1;
  canname+=canname+"-ext";
  TCanvas *c1f=new TCanvas(canname.c_str(),"CANFITFUNC", 900,900);
  c1f->SetFillColor(kWhite);
  c1f->SetBorderMode(0);
  c1f->cd();
  f1->GetXaxis()->SetTitle("M_{H}");
  f1->GetYaxis()->SetTitle("#varepsilon");
  f1->GetYaxis()->SetTitleOffset(1.0);
  f1->Draw("L");
  gr_eff->Draw("Psame");
  c1f->SaveAs( ("eff_param_newHLT_"+btag_str +chan_str+"-extrange.root").c_str()  );
  c1f->SaveAs( ("eff_param_newHLT_"+btag_str +chan_str+"-extrange.C").c_str()  );
  // gr_eff->Write();

  TFile *outf1=new TFile("./convertedTree_LP_20110811.root","UPDATE");
  outf1->cd();
  gr_eff->Write();
  f1->Write();
  outf1->Close();
  delete gr_eff;
  delete c1f;
  delete f1;
  return fit_res;
}//end eff_fit


vector<double> calculate_CBpars(int ibtag,double mH){
  vector<double>outpars;
  double CB_mean,CB_sigma,CB_alpha,CB_n,poly0, poly1;
  if(ibtag==0){

    CB_mean = -2.92969 + 0.0444168*mH;
    CB_sigma = -1.17027 + .0394395*mH;
    CB_alpha = 1.6114;
    CB_n = .618079 + .000493633*mH;
    poly0 = -2.08918 - .158074*mH + .000334758*mH*mH;
    poly1 = 0.11238;
  }
  else if(ibtag==1){
    CB_mean= -4.63689+0.0514951*mH;
    CB_sigma= -8.64727+0.066107*mH;
    CB_alpha= .477979+0.00485747*mH- 0.00000542886*mH*mH;
    CB_n= 0.264879+.00105916*mH;
    poly0= -6.64749-0.13105*mH+0.000297844*mH*mH;
    poly1= 0.110152;
  }
  else if(ibtag==2){ 

    CB_mean = -2.77113 + .0442444*mH;
    CB_sigma = -1.31594 + .0375571*mH ;
    CB_alpha = 1.07316 + .0014977*mH;
    CB_n = 2.57449 - .0029208*mH;
    poly0 = -21.0;
    poly1 = .1095;
  }
  else cout<<"How many btags do you think to have??? "<<ibtag<<endl;

  outpars.push_back(CB_mean);
  outpars.push_back(CB_sigma);
  outpars.push_back(CB_alpha);
  outpars.push_back(CB_n);

  outpars.push_back(poly0);
  outpars.push_back(poly1);

  return outpars;
}//end vector<double> calculate_CBpars


void extract_exp_sig_yields(float mymass,float mywidth){
  int im=get_mass_index(mymass);
  if(im<0)return;

 ///////////CHANGE THIS DIR NAME ALSO IN get_gen_yields() 
  //string myDir="/afs/cern.ch/user/s/sbologne/scratch0/CMSSW/CMSSW_4_2_4/src/HiggsAnalysis/CombinedLimit/test/rotatedEPSForLP/treeFromFrancesco/";
  string myDir="/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/HZZlljj/HZZlljjAnalyzer/test/analysis/FIT/";
  //1nvfb_trigfix//";//dir with Francesco's tree //960invpb/

  //signal file
  int powheg_jhu = get_signal_gen(mymass);
  string prefix;
  string suffix;
  if(powheg_jhu==0){
    prefix="HZZlljjRM_GluGluToHToZZTo2L2Q_M-";
    suffix="_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root";
  }
  else{
    prefix="HZZlljjRM_SMHiggsToZZTo2L2Q_M-";
    suffix="_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2_optLD_looseBTags_v2_ALL.root";
  }
  std::ostringstream ossm;
  ossm<<mass[im];
  string filename=myDir+prefix+ossm.str()+suffix;
  cout << "Opening file: " << filename << " for calculating number of selected signal events" << endl;
  TFile *f=new TFile(filename.c_str(),"READ");
  TTree *t=(TTree*)f->Get("tree_passedEvents");

  //bkg file
  /*string filebkg=myDir+"HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root";
  TFile *fbkg=new TFile(filebkg.c_str(),"READ");
  TTree *tbkg=(TTree*)f->Get("tree_passedEvents");*/

  string baseline_sel="mZjj>75&&mZjj<105&&mZZ>183.0";
  // cout<<"MyWidth "<<mymass<<" OLD: "<<mywidth<<flush;
  //mywidth = myCSW->HiggsWidth(0,mymass);
  // cout<<"  NEW: "<<mywidth<<endl;
  double effWidth=sqrt(mywidth*mywidth+100);  // effective width used for defining window
  cout << "MYWIDTH: " << mywidth << endl;
  cout << "EFFECTIVE WIDTH: " << effWidth << endl;
  double lowMZZ=99.;
  double hiMZZ=101.0; 

  if((mymass-10*effWidth)<183.0)  lowMZZ=183.0;
  else lowMZZ=mymass-10*effWidth;

  if((mymass+10*effWidth)>800.0)  hiMZZ=800.;
  else hiMZZ=mymass+10*effWidth;


  //float lowMZZ= mymass-sqrt(mywidth*mywidth+100.0)*10.0;
  //float hiMZZ = mymass+sqrt(mywidth*mywidth+100.0)*10.0;
  //if(lowMZZ<183.0)lowMZZ=183.0;
  //if(hiMZZ>800.0)hiMZZ=800.0;
  std::ostringstream oss1;
  oss1<<lowMZZ;
  std::ostringstream oss2;
  oss2<<hiMZZ;
  string massCut="&&mZZ>"+oss1.str()+" && mZZ<"+oss2.str();
  baseline_sel+=massCut;

 
  for(int ib=0;ib<nbtag;ib++){
    string sel="dumm";
    if(ib==0)sel=baseline_sel+"&&nBTags==0";
    else if(ib==1)sel=baseline_sel+"&&nBTags==1";
    else if(ib==2)sel=baseline_sel+"&&nBTags==2";
    string sel2="dumsel2";
    for(int ich=0;ich<nchan;ich++){
      cout << "CUTSTRING - " << sel << endl;
      if(ich==0)sel2="("+sel+"&&leptType==1.0)*HLTSF";//*HLTSF  electrons, trees have different convention than this code    
      else sel2="("+sel+"&&leptType==0.0)*HLTSF";//*HLTSF
      //      cout<<"M"<<mymass<<", "<<ib<<"b, "<<(ich==0 ? "ee" : "mm") <<"\tSelection string "<<sel2.c_str()<<endl;


      float nexp=float(t->Draw("mZZ>>mzzrescaled(400,183.0,1200.0)",sel2.c_str(),"goff"));//sel2 contains tighter mzz cut
      TH1F *mzztmp=(TH1F*)gDirectory->Get("mzzrescaled");
      nexp=float( mzztmp->Integral()  );
      exp_sig_yields[ib][ich][im]=nexp;
      cout << "number of selected events for " << mymass << " is: " << nexp << endl;
      delete mzztmp;
      //cout<<"MH-"<<mymass<<" "<<ib<<"b, ee only -> "<<im<<" Selecting : "<<sel2.c_str()<<"  Exp Sig "<<exp_sig_yields[ib][ich][im]<<"  ExpBkg "<<exp_bkg[ib][ich]<<endl;
    }
  }//end loop on btag

  delete t;delete f;


}//extract_exp_sig_yields

void get_gen_yields(){
  //string myDir="/afs/cern.ch/user/s/sbologne/scratch0/CMSSW/CMSSW_4_2_4/src/HiggsAnalysis/CombinedLimit/test/rotatedEPSForLP/treeFromFrancesco/";
  string myDir="/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/HZZlljj/HZZlljjAnalyzer/test/analysis/FIT/";
  //signal file
  for(int im=0;im<nmass;im++){
   int powheg_jhu = get_signal_gen(mass[im]);
  string prefix;
  string suffix;
  if(powheg_jhu==0){
    prefix="HZZlljjRM_GluGluToHToZZTo2L2Q_M-";
    suffix="_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root";
  }
  else{
    prefix="HZZlljjRM_SMHiggsToZZTo2L2Q_M-";
    suffix="_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2_optLD_looseBTags_v2_ALL.root";
  }
   std::ostringstream ossm;
    ossm<<mass[im];
    string filename=myDir+prefix+ossm.str()+suffix;

    //-------------------- AJW ---------------
    cout << "Opening: " << filename << " for extracting number of generated events." << endl;
    TFile *f=new TFile(filename.c_str(),"READ");
    TH1F *h=(TH1F*)f->Get("nCounter");
    gen_sig_yields[im]=h->GetBinContent(1);
    cout << "the generated number of events for " << mass[im] << " is: " << gen_sig_yields[im] << endl;
    if(gen_sig_yields[im]<=0.0){
      cout<<"WARNING!!! # gen events for mass "<<mass[im]<<"  == "<<gen_sig_yields[im]<<"  . Setting to -1."<<endl;
      gen_sig_yields[im]=-1.0;
    }
  }

}

int get_signal_gen(float mymass){
  int genID=-1;//0=powheg; 1=JHUGEN
  /*
 if(mymass==190.0||mymass==200.0||mymass==210.0||
    mymass==230.0||mymass==300.0||mymass==425.0||
    mymass==525.0||mymass==550.0||mymass==575.0||
    mymass==600.0) genID=0;
 else genID=1;
  */

  if(mymass==250.0||mymass==300.0||mymass==350.0||mymass==400.0||
     mymass==450.0||mymass==500.0) genID=1;
  else genID=0;
  return genID;

}
int get_mass_index(float mymass){
  int index=-1;
  //  cout<<"Mass "<<mymass<<flush;
  for(int im=0;im<nmass;im++){
    if(mymass==mass[im]){
      index=im;
      break;
    }
  }
  //  if(index>=0)cout<<" found"<<endl;
  // else cout<<" NOT found"<<endl;
  return index;

}
