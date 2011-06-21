# --------------------------------------------------
# USAGE: source finalize_all/sh [selType] [leptType]
# leptType can be "ELE", "MU" or "ALL"
# --------------------------------------------------

if( $2 == "ELE" ) then
./finalize_HZZlljj DoubleMu_Run2011A $1 $2
else if( $2 == "MU" ) then
./finalize_HZZlljj DoubleElectron_Run2011A $1 $2
else
./finalize_HZZlljj DATA_Run2011A $1 $2
endif
./finalize_HZZlljj SMHiggsToZZTo2L2Q_M-250_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2 $1 $2
./finalize_HZZlljj SMHiggsToZZTo2L2Q_M-300_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2 $1 $2
./finalize_HZZlljj SMHiggsToZZTo2L2Q_M-350_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2 $1 $2
./finalize_HZZlljj SMHiggsToZZTo2L2Q_M-400_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2 $1 $2
./finalize_HZZlljj SMHiggsToZZTo2L2Q_M-450_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2 $1 $2
./finalize_HZZlljj SMHiggsToZZTo2L2Q_M-500_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2 $1 $2
./finalize_HZZlljj GluGluToHToZZTo2L2Q_M-600_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2 $1 $2
./finalize_HZZlljj VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11_v2 $1 $2
./finalize_HZZlljj TT_TW_TuneZ2_7TeV-pythia6-tauola_Spring11_v2 $1 $2
./finalize_HZZlljj ZJets_alpgen_TuneZ2_Spring11_v2 $1 $2
./finalize_HZZlljj ZCC_alpgen_TuneZ2_Spring11_v2 $1 $2
./finalize_HZZlljj ZBB_alpgen_TuneZ2_Spring11_v2 $1 $2
