# --------------------------------------------------
# USAGE: source finalize_all/sh [selType] [leptType]
# leptType can be "ELE", "MU" or "ALL"
# --------------------------------------------------

if( $2 == "ELE" ) then
./finalize_HZZlljj Electron_Nov4ReReco_PU $1 $2
else if( $2 == "MU" ) then
./finalize_HZZlljj Mu_Nov4ReReco_PU $1 $2
else
./finalize_HZZlljj EleMu_Nov4ReReco_PU $1 $2
endif
./finalize_HZZlljj SMHiggsToZZTo2L2Q_M-250_7TeV-jhu-pythia6 $1 $2
./finalize_HZZlljj SMHiggsToZZTo2L2Q_M-300_7TeV-jhu-pythia6 $1 $2
./finalize_HZZlljj SMHiggsToZZTo2L2Q_M-350_7TeV-jhu-pythia6 $1 $2
./finalize_HZZlljj SMHiggsToZZTo2L2Q_M-400_7TeV-jhu-pythia6 $1 $2
./finalize_HZZlljj SMHiggsToZZTo2L2Q_M-450_7TeV-jhu-pythia6 $1 $2
./finalize_HZZlljj SMHiggsToZZTo2L2Q_M-500_7TeV-jhu-pythia6 $1 $2
./finalize_HZZlljj VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10 $1 $2
./finalize_HZZlljj TTJets_TuneZ2_7TeV-madgraph-tauola_Fall10 $1 $2
#./finalize_HZZlljj ZJets_madgraph $1
./finalize_HZZlljj ZJets_alpgen_TuneZ2_Fall10 $1 $2
