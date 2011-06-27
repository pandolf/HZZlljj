# --------------------------------------------------
# USAGE: source finalize_all/sh [selType] [leptType]
# leptType can be "ELE", "MU" or "ALL"
# --------------------------------------------------

#if( $2 == "ELE" ) then
#./finalize_HZZlljjRM DoubleMu_Run2011A $1 $2
#else if( $2 == "MU" ) then
#./finalize_HZZlljjRM DoubleElectron_Run2011A $1 $2
#else
#./finalize_HZZlljjRM DATA_Run2011A_v2_Sub2 $1 $2
#endif
./finalize_HZZlljjRM GluGluToHToZZTo2L2Q_M-190_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1 $1 $2
./finalize_HZZlljjRM GluGluToHToZZTo2L2Q_M-200_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1 $1 $2
./finalize_HZZlljjRM GluGluToHToZZTo2L2Q_M-210_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1 $1 $2
./finalize_HZZlljjRM GluGluToHToZZTo2L2Q_M-230_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1 $1 $2
#./finalize_HZZlljjRM SMHiggsToZZTo2L2Q_M-300_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2 $1 $2
#./finalize_HZZlljjRM SMHiggsToZZTo2L2Q_M-350_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2 $1 $2
#./finalize_HZZlljjRM SMHiggsToZZTo2L2Q_M-400_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2 $1 $2
#./finalize_HZZlljjRM SMHiggsToZZTo2L2Q_M-450_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2 $1 $2
#./finalize_HZZlljjRM SMHiggsToZZTo2L2Q_M-500_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2 $1 $2
#./finalize_HZZlljjRM GluGluToHToZZTo2L2Q_M-600_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2 $1 $2
#./finalize_HZZlljjRM VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11_v2 $1 $2
#./finalize_HZZlljjRM TT_TW_TuneZ2_7TeV-pythia6-tauola_Spring11_v2 $1 $2
#./finalize_HZZlljjRM ZJets_alpgen_TuneZ2_Spring11_v2 $1 $2
#./finalize_HZZlljjRM ZCC_alpgen_TuneZ2_Spring11_v2 $1 $2
#./finalize_HZZlljjRM ZBB_alpgen_TuneZ2_Spring11_v2 $1 $2
