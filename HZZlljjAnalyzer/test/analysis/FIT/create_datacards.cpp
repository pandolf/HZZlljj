#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>




void create_singleDatacard( float mass, float lumi, const std::string& leptType_str, int nbtags );


int main( int argc, char* argv[] ) {

  if( argc!=2 ) {
    std::cout << "USAGE: ./create_datacards [dataset]" << std::endl;
    exit(113);
  }

  
  std::string dataset(argv[1]);

  float lumi_ELE;
  float lumi_MU;
  if( dataset=="LP11") {
    lumi_ELE=1556.; //pb^-1
    lumi_MU =1615.; //pb^-1
  } else if( dataset=="Run2011A_FULL" ) {
    lumi_ELE=2100.; //pb^-1
    lumi_MU =2100.; //pb^-1
  } else {
    std::cout << "Unknown dataset '" << dataset << "'. Exiting." << std::endl;
    exit(333);
  }


  std::ifstream ifs("masses.txt");
  
  while( ifs.good() ) {
    
    float mass;
    ifs >> mass;

    std::cout << "+++ MASS: " << mass << std::endl;

    char mkdir_command[100];
    sprintf( mkdir_command, "mkdir -p datacardsPROVA/%.0f", mass);
    system(mkdir_command);

    create_singleDatacard( mass, lumi_ELE, "ELE", 0);
    create_singleDatacard( mass, lumi_ELE, "ELE", 1);
    create_singleDatacard( mass, lumi_ELE, "ELE", 2);
    create_singleDatacard( mass, lumi_MU,  "MU", 0);
    create_singleDatacard( mass, lumi_MU,  "MU", 1);
    create_singleDatacard( mass, lumi_MU,  "MU", 2);

  } //while masses

  return 0;

}



void create_singleDatacard( float mass, float lumi, const std::string& leptType_str, int nbtags ) {

  if( leptType_str!="ELE" && leptType_str!="MU" ) {
    std::cout << "Unkown Lept Type '" << leptType_str << "'. Exiting." << std::endl;
    exit(12333);
  }

  
  std::string leptType_forDatacard = (leptType_str=="ELE") ? "ee" : "mm";
  char suffix[100];
  sprintf( suffix, "%s%db", leptType_forDatacard.c_str(), nbtags);
  std::string suffix_str(suffix);

  char datacardName[400];
  sprintf( datacardName, "datacardsPROVA/%.0f/hzz2l2q_%s.%.0f.txt", mass, suffix, mass);


  std::ofstream ofs(datacardName);
  ofs << "# Simple counting experiment, with one signal and one background process" << std::endl;
  ofs << "#imax 1  number of channels" << std::endl;
  ofs << "#jmax 1  number of backgrounds" << std::endl;
  ofs << "#kmax *  number of nuisance parameters (sources of systematical uncertainties)" << std::endl;
  ofs << "------------ " << std::endl;
  ofs << "shapes ggH CMS_hzz2l2q_" << suffix_str << " hzz2l2q_" << suffix_str << ".input.root  w:signal" << std::endl;
  ofs << "shapes VBF CMS_hzz2l2q_" << suffix_str << " hzz2l2q_" << suffix_str << ".input.root  w:signal " << std::endl;
  ofs << "shapes background CMS_hzz2l2q_" << suffix_str << " hzz2l2q_" << suffix_str << ".input.root w:background " << std::endl;
  ofs << "shapes data_obs   CMS_hzz2l2q_" << suffix_str << " hzz2l2q_" << suffix_str << ".input.root w:dataset_obs" << std::endl;
  ofs << "------------ " << std::endl;
  ofs << "# we have just one channel, in which we observe 0 events " << std::endl;
  ofs << "bin         CMS_hzz2l2q_" << suffix << std::endl;






  ofs.close();

}
