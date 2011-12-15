CC = g++
CFLAGS = -Wall -c -g

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/

ROOTSYS  ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs` 

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit


merge_and_setWeights: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/merge_and_setWeights.cpp
	$(CC) -Wall $(INCLUDES) -o merge_and_setWeights $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/merge_and_setWeights.cpp $(ROOTFLAG) $(EXTRALIBS)

merge_and_setWeights_HWWlvjj: merge_and_setWeights_HWWlvjj.cpp
	$(CC) -Wall -o merge_and_setWeights_HWWlvjj merge_and_setWeights_HWWlvjj.cpp $(ROOTFLAG)

#do2ndLevel_HZZlljj: Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o QGLikelihoodCalculator.o do2ndLevel_HZZlljj.cpp
#	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_HZZlljj do2ndLevel_HZZlljj.cpp Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o QGLikelihoodCalculator.o $(ROOTFLAG) -lCore -lTMVA

do2ndLevel_HZZlljj: Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o AnalysisElectron.o AnalysisMuon.o PUWeight.o do2ndLevel_HZZlljj.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_HZZlljj do2ndLevel_HZZlljj.cpp Ntp1Analyzer.o PUWeight.o Ntp1Analyzer_HZZlljj.o AnalysisElectron.o AnalysisMuon.o $(ROOTFLAG) 

do2ndLevel_DumpNvertex: Ntp1Analyzer.o Ntp1Analyzer_DumpNvertex.o AnalysisElectron.o AnalysisMuon.o PUWeight.o do2ndLevel_DumpNvertex.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_DumpNvertex do2ndLevel_DumpNvertex.cpp Ntp1Analyzer.o PUWeight.o Ntp1Analyzer_DumpNvertex.o AnalysisElectron.o AnalysisMuon.o $(ROOTFLAG) 

do2ndLevel_HWWlvjj: Ntp1Analyzer.o Ntp1Analyzer_HWWlvjj.o AnalysisElectron.o AnalysisMuon.o Utils.o do2ndLevel_HWWlvjj.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_HWWlvjj do2ndLevel_HWWlvjj.cpp Ntp1Analyzer.o Ntp1Analyzer_HWWlvjj.o AnalysisElectron.o AnalysisMuon.o Utils.o $(ROOTFLAG)

do2ndLevel_WZ: Ntp1Analyzer.o Ntp1Analyzer_WZ.o do2ndLevel_WZ.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_WZ do2ndLevel_WZ.cpp Ntp1Analyzer.o Ntp1Analyzer_WZ.o $(ROOTFLAG) 

do2ndLevel_WZ_DATA: Ntp1Analyzer.o Ntp1Analyzer_WZ.o do2ndLevel_WZ_DATA.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_WZ_DATA do2ndLevel_WZ_DATA.cpp Ntp1Analyzer.o Ntp1Analyzer_WZ.o $(ROOTFLAG) 

do2ndLevel_HZZlljj_DATA: Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o AnalysisElectron.o AnalysisMuon.o PUWeight.o do2ndLevel_HZZlljj_DATA.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_HZZlljj_DATA do2ndLevel_HZZlljj_DATA.cpp Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o AnalysisElectron.o AnalysisMuon.o PUWeight.o $(ROOTFLAG)

do2ndLevel_LeptIsoStudies: Ntp1Analyzer.o Ntp1Analyzer_LeptIsoStudies.o AnalysisElectron.o AnalysisMuon.o do2ndLevel_LeptIsoStudies.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_LeptIsoStudies do2ndLevel_LeptIsoStudies.cpp Ntp1Analyzer.o Ntp1Analyzer_LeptIsoStudies.o AnalysisElectron.o AnalysisMuon.o $(ROOTFLAG) 



do2ndLevel_HWWlvjj_DATA: Ntp1Analyzer.o Ntp1Analyzer_HWWlvjj.o AnalysisElectron.o AnalysisMuon.o Utils.o do2ndLevel_HWWlvjj_DATA.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_HWWlvjj_DATA do2ndLevel_HWWlvjj_DATA.cpp Ntp1Analyzer.o Ntp1Analyzer_HWWlvjj.o AnalysisElectron.o AnalysisMuon.o Utils.o $(ROOTFLAG)

finalize_HWWlvjj: Ntp1Finalizer.o Ntp1Finalizer_HWWlvjj.o finalize_HWWlvjj.cpp fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitConstraintMGaus.o TFitConstraintEp.o TFitParticleEtEtaPhi.o TKinFitter.o QGLikelihoodCalculator.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o DiJetKinFitter.o LeptonNeutrinoKinFitter.o GlobalFitter.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood -o finalize_HWWlvjj finalize_HWWlvjj.cpp Ntp1Finalizer.o Ntp1Finalizer_HWWlvjj.o fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitConstraintMGaus.o TFitConstraintEp.o TFitParticleEtEtaPhi.o TKinFitter.o QGLikelihoodCalculator.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o DiJetKinFitter.o LeptonNeutrinoKinFitter.o GlobalFitter.o $(ROOTFLAG) -lMinuit


finalize_HZZlljj: Ntp1Finalizer.o Ntp1Finalizer_HZZlljj.o finalize_HZZlljj.cpp fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitParticleEtEtaPhi.o TKinFitter.o QGLikelihoodCalculator.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o DiJetKinFitter.o AnalysisJet.o BTagSFUtil.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood -o finalize_HZZlljj finalize_HZZlljj.cpp Ntp1Finalizer.o Ntp1Finalizer_HZZlljj.o fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitParticleEtEtaPhi.o TKinFitter.o QGLikelihoodCalculator.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o DiJetKinFitter.o AnalysisJet.o BTagSFUtil.o $(ROOTFLAG) -lMinuit

finalize_HZZlljjRM: Ntp1Finalizer.o Ntp1Finalizer_HZZlljjRM.o finalize_HZZlljjRM.cpp fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitParticleEtEtaPhi.o TKinFitter.o QGLikelihoodCalculator.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o DiJetKinFitter.o AnalysisJet.o BTagSFUtil.o PUWeight.o
	$(CC) -Wall $(INCLUDES) -o finalize_HZZlljjRM finalize_HZZlljjRM.cpp Ntp1Finalizer.o Ntp1Finalizer_HZZlljjRM.o fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitParticleEtEtaPhi.o TKinFitter.o QGLikelihoodCalculator.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o DiJetKinFitter.o AnalysisJet.o BTagSFUtil.o PUWeight.o $(ROOTFLAG) $(EXTRALIBS)

finalize_WZ: Ntp1Finalizer.o Ntp1Finalizer_WZ.o finalize_WZ.cpp fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitParticleEtEtaPhi.o TKinFitter.o QGLikelihoodCalculator.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o DiJetKinFitter.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood -o finalize_WZ finalize_WZ.cpp Ntp1Finalizer.o Ntp1Finalizer_WZ.o fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitParticleEtEtaPhi.o TKinFitter.o QGLikelihoodCalculator.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o DiJetKinFitter.o $(ROOTFLAG) -lMinuit

finalize_ZGamma: Ntp1Finalizer.o Ntp1Finalizer_ZGamma.o finalize_ZGamma.cpp fitTools.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o finalize_ZGamma finalize_ZGamma.cpp Ntp1Finalizer.o Ntp1Finalizer_ZGamma.o fitTools.o $(ROOTFLAG)

finalize_QG: Ntp1Finalizer.o Ntp1Finalizer_QG.o finalize_QG.cpp fitTools.o AnalysisJet.o BTagSFUtil.o
	$(CC) -Wall $(INCLUDES) -o finalize_QG finalize_QG.cpp Ntp1Finalizer.o Ntp1Finalizer_QG.o fitTools.o AnalysisJet.o BTagSFUtil.o $(ROOTFLAG) $(EXTRALIBS)

finalize_QGTMVA: Ntp1Finalizer.o Ntp1Finalizer_QGTMVA.o finalize_QGTMVA.cpp fitTools.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o finalize_QGTMVA finalize_QGTMVA.cpp Ntp1Finalizer.o Ntp1Finalizer_QGTMVA.o fitTools.o $(ROOTFLAG) -lMinuit

finalize_ComputeQGLikelihood: Ntp1Finalizer.o Ntp1Finalizer_ComputeQGLikelihood.o QGLikelihoodCalculator.o finalize_ComputeQGLikelihood.cpp fitTools.o AnalysisJet.o BTagSFUtil.o
	$(CC) -Wall $(INCLUDES) -o finalize_ComputeQGLikelihood finalize_ComputeQGLikelihood.cpp Ntp1Finalizer.o Ntp1Finalizer_ComputeQGLikelihood.o QGLikelihoodCalculator.o fitTools.o AnalysisJet.o BTagSFUtil.o $(ROOTFLAG) $(EXTRALIBS)

finalize_GG: Ntp1Finalizer.o Ntp1Finalizer_GG.o finalize_GG.cpp fitTools.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o finalize_GG finalize_GG.cpp Ntp1Finalizer.o Ntp1Finalizer_GG.o fitTools.o $(ROOTFLAG) -lMinuit

finalize_KinFit: Ntp1Finalizer.o Ntp1Finalizer_KinFit.o finalize_KinFit.cpp fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitConstraintMGaus.o TFitConstraintEp.o TFitParticleEtEtaPhi.o TKinFitter.o QGLikelihoodCalculator.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood -o finalize_KinFit finalize_KinFit.cpp Ntp1Finalizer.o Ntp1Finalizer_KinFit.o fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitConstraintMGaus.o TFitConstraintEp.o TFitParticleEtEtaPhi.o TKinFitter.o QGLikelihoodCalculator.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o $(ROOTFLAG) -lMinuit


do2ndLevel_JetStudies: Ntp1Analyzer.o Ntp1Analyzer_JetStudies.o 
	$(CC) -Wall -o do2ndLevel_JetStudies do2ndLevel_JetStudies.cpp Ntp1Analyzer.o Ntp1Analyzer_JetStudies.o $(ROOTFLAG)

do2ndLevel_Zbb: Ntp1Analyzer.o Ntp1Analyzer_Zbb.o do2ndLevel_Zbb.cpp
	$(CC) -Wall -o do2ndLevel_Zbb do2ndLevel_Zbb.cpp Ntp1Analyzer.o Ntp1Analyzer_Zbb.o $(ROOTFLAG)

do2ndLevel_ZGamma: Ntp1Analyzer.o Ntp1Analyzer_ZGamma.o 
	$(CC) -Wall -o do2ndLevel_ZGamma do2ndLevel_ZGamma.cpp Ntp1Analyzer.o Ntp1Analyzer_ZGamma.o $(ROOTFLAG)

do2ndLevel_TMVA: Ntp1Analyzer.o Ntp1Analyzer_TMVA.o do2ndLevel_TMVA.cpp QGLikelihoodCalculator.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitParticleEtEtaPhi.o TKinFitter.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o
	$(CC) -Wall -o do2ndLevel_TMVA do2ndLevel_TMVA.cpp Ntp1Analyzer.o Ntp1Analyzer_TMVA.o QGLikelihoodCalculator.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitParticleEtEtaPhi.o TKinFitter.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o $(ROOTFLAG) 

do2ndLevel_QG: Ntp1Analyzer.o Ntp1Analyzer_QG.o do2ndLevel_QG.cpp
	$(CC) -Wall -o do2ndLevel_QG do2ndLevel_QG.cpp Ntp1Analyzer.o Ntp1Analyzer_QG.o $(ROOTFLAG)

do2ndLevel_DiJet: Ntp1Analyzer.o Ntp1Analyzer_DiJet.o do2ndLevel_DiJet.cpp
	$(CC) -Wall -o do2ndLevel_DiJet do2ndLevel_DiJet.cpp Ntp1Analyzer.o Ntp1Analyzer_DiJet.o $(ROOTFLAG)

drawHWWlvjj: DrawBase.o fitTools.o drawHWWlvjj.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawHWWlvjj drawHWWlvjj.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawQGCompareDATA: DrawBase.o fitTools.o drawQGCompareDATA.cpp
	$(CC) -Wall $(INCLUDES) -o drawQGCompareDATA drawQGCompareDATA.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

drawHZZlljjRM: DrawBase.o fitTools.o drawHZZlljjRM.cpp
	$(CC) -Wall $(INCLUDES) -o drawHZZlljjRM drawHZZlljjRM.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

drawHZZlljjRM_DATA: DrawBase.o fitTools.o drawHZZlljjRM_DATA.cpp
	$(CC) -Wall $(INCLUDES) -o drawHZZlljjRM_DATA drawHZZlljjRM_DATA.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

drawHZZlljjRM_DATA_Check1TeV: DrawBase.o fitTools.o drawHZZlljjRM_DATA_Check1TeV.cpp
	$(CC) -Wall $(INCLUDES) -o drawHZZlljjRM_DATA_Check1TeV drawHZZlljjRM_DATA_Check1TeV.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

drawHZZlljjRM_JHU_vs_Powheg: DrawBase.o fitTools.o drawHZZlljjRM_JHU_vs_Powheg.cpp
	$(CC) -Wall $(INCLUDES) -o drawHZZlljjRM_JHU_vs_Powheg drawHZZlljjRM_JHU_vs_Powheg.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

drawHZZlljj_kinfit: DrawBase.o fitTools.o drawHZZlljj_kinfit.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawHZZlljj_kinfit drawHZZlljj_kinfit.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawHZZlljj_DATA: DrawBase.o fitTools.o drawHZZlljj_DATA.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawHZZlljj_DATA drawHZZlljj_DATA.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawKinFit: DrawBase.o fitTools.o drawKinFit.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf -o drawKinFit drawKinFit.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawWZ_DATA: DrawBase.o fitTools.o drawWZ_DATA.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf -o drawWZ_DATA drawWZ_DATA.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawMCHistos: DrawBase.o fitTools.o drawMCHistos.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawMCHistos drawMCHistos.cpp DrawBase.o fitTools.o $(ROOTFLAG)

drawJetStudies: DrawBase.o fitTools.o drawJetStudies.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawJetStudies drawJetStudies.cpp DrawBase.o fitTools.o $(ROOTFLAG)

drawZGamma: DrawBase.o fitTools.o drawZGamma.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawZGamma drawZGamma.cpp DrawBase.o fitTools.o $(ROOTFLAG)

drawZGamma_DATA: DrawBase.o fitTools.o drawZGamma_DATA.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawZGamma_DATA drawZGamma_DATA.cpp DrawBase.o fitTools.o $(ROOTFLAG)

drawKinFitBias: DrawBase.o fitTools.o drawKinFitBias.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawKinFitBias drawKinFitBias.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawQG: DrawBase.o fitTools.o drawQG.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawQG drawQG.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

drawRej_vs_Eff: DrawBase.o fitTools.o drawRej_vs_Eff.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawRej_vs_Eff drawRej_vs_Eff.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawRej_vs_Eff_LD: DrawBase.o fitTools.o drawRej_vs_Eff_LD.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawRej_vs_Eff_LD drawRej_vs_Eff_LD.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawJetStudies_H: DrawBase.o fitTools.o drawJetStudies_H.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawJetStudies_H drawJetStudies_H.cpp DrawBase.o fitTools.o $(ROOTFLAG)

draw_Fall10_vs_Winter10: DrawBase.o fitTools.o draw_Fall10_vs_Winter10.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o draw_Fall10_vs_Winter10 draw_Fall10_vs_Winter10.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawThresholds_vs_Mass: DrawBase.o fitTools.o drawThresholds_vs_Mass.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawThresholds_vs_Mass drawThresholds_vs_Mass.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawThresholdsUL_vs_Mass: DrawBase.o fitTools.o drawThresholdsUL_vs_Mass.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawThresholdsUL_vs_Mass drawThresholdsUL_vs_Mass.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawFall10_vs_Spring11: DrawBase.o fitTools.o drawFall10_vs_Spring11.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawFall10_vs_Spring11 drawFall10_vs_Spring11.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawLeptIsoStudies: DrawBase.o fitTools.o drawLeptIsoStudies.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawLeptIsoStudies drawLeptIsoStudies.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawAlpha: DrawBase.o fitTools.o drawAlpha.cpp
	$(CC) -Wall $(INCLUDES) -o drawAlpha drawAlpha.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)


doLatexYieldTable: doLatexYieldTable.cpp
	$(CC) -Wall -o doLatexYieldTable doLatexYieldTable.cpp $(ROOTFLAG)

fitSidebands: fitSidebands.cpp DrawBase.o fitTools.o RooFermi.o RooCruijffPdf.o
	$(CC) -Wall $(INCLUDES) -o fitSidebands fitSidebands.cpp DrawBase.o fitTools.o RooFermi.o RooCruijffPdf.o $(ROOTFLAG) $(EXTRALIBS)



Ntp1Analyzer.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Analyzer.C
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Analyzer.C $(ROOTFLAG)

Ntp1Analyzer_HZZlljj.o: Ntp1Analyzer_HZZlljj.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/emanuele/CommonTools Ntp1Analyzer_HZZlljj.C $(ROOTFLAG)

Ntp1Analyzer_DumpNvertex.o: Ntp1Analyzer_DumpNvertex.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/emanuele/CommonTools Ntp1Analyzer_DumpNvertex.C $(ROOTFLAG)

Ntp1Analyzer_HWWlvjj.o: Ntp1Analyzer_HWWlvjj.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Analyzer_HWWlvjj.C $(ROOTFLAG)

Ntp1Analyzer_WZ.o: Ntp1Analyzer_WZ.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Analyzer_WZ.C $(ROOTFLAG)

Ntp1Analyzer_JetStudies.o: Ntp1Analyzer_JetStudies.C
	$(CC) $(CFLAGS) Ntp1Analyzer_JetStudies.C $(ROOTFLAG)

Ntp1Analyzer_ZGamma.o: Ntp1Analyzer_ZGamma.C
	$(CC) $(CFLAGS) Ntp1Analyzer_ZGamma.C $(ROOTFLAG)

Ntp1Analyzer_TMVA.o: Ntp1Analyzer_TMVA.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/ Ntp1Analyzer_TMVA.C $(ROOTFLAG)

Ntp1Analyzer_QG.o: Ntp1Analyzer_QG.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Analyzer_QG.C $(ROOTFLAG)

Ntp1Analyzer_DiJet.o: Ntp1Analyzer_DiJet.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Analyzer_DiJet.C $(ROOTFLAG)

Ntp1Analyzer_Zbb.o: Ntp1Analyzer_Zbb.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Analyzer_Zbb.C $(ROOTFLAG)

Ntp1Finalizer.o: Ntp1Finalizer.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/ Ntp1Finalizer.C $(ROOTFLAG)

Ntp1Finalizer_HZZlljj.o: Ntp1Finalizer_HZZlljj.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood Ntp1Finalizer_HZZlljj.C $(ROOTFLAG)

Ntp1Finalizer_HZZlljjRM.o: Ntp1Finalizer_HZZlljjRM.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood Ntp1Finalizer_HZZlljjRM.C $(ROOTFLAG)

Ntp1Finalizer_HWWlvjj.o: Ntp1Finalizer_HWWlvjj.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood Ntp1Finalizer_HWWlvjj.C $(ROOTFLAG)

Ntp1Finalizer_WZ.o: Ntp1Finalizer_WZ.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood Ntp1Finalizer_WZ.C $(ROOTFLAG)

Ntp1Finalizer_ZGamma.o: Ntp1Finalizer_ZGamma.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/  Ntp1Finalizer_ZGamma.C $(ROOTFLAG)

Ntp1Finalizer_QG.o: Ntp1Finalizer_QG.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/  Ntp1Finalizer_QG.C $(ROOTFLAG)

Ntp1Finalizer_QGTMVA.o: Ntp1Finalizer_QGTMVA.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/  Ntp1Finalizer_QGTMVA.C $(ROOTFLAG)

Ntp1Finalizer_ComputeQGLikelihood.o: Ntp1Finalizer_ComputeQGLikelihood.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf Ntp1Finalizer_ComputeQGLikelihood.C $(ROOTFLAG)

Ntp1Finalizer_GG.o: Ntp1Finalizer_GG.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/  Ntp1Finalizer_GG.C $(ROOTFLAG)

Ntp1Finalizer_KinFit.o: Ntp1Finalizer_KinFit.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/ Ntp1Finalizer_KinFit.C $(ROOTFLAG)

TFitParticleEtEtaPhi.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitParticleEtEtaPhi.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitParticleEtEtaPhi.cc $(ROOTFLAG)

TAbsFitParticle.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TAbsFitParticle.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TAbsFitParticle.cc $(ROOTFLAG)

TFitConstraintM.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintM.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintM.cc $(ROOTFLAG)

TFitConstraintMGaus.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintMGaus.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintMGaus.cc $(ROOTFLAG)

TFitConstraintEp.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintEp.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintEp.cc $(ROOTFLAG)

TAbsFitConstraint.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TAbsFitConstraint.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TAbsFitConstraint.cc $(ROOTFLAG)

TKinFitter.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TKinFitter.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TKinFitter.cc $(ROOTFLAG)

DiJetKinFitter.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/DiJetKinFitter.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/DiJetKinFitter.cc $(ROOTFLAG)

LeptonNeutrinoKinFitter.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/LeptonNeutrinoKinFitter.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/LeptonNeutrinoKinFitter.cc $(ROOTFLAG)

GlobalFitter.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/GlobalFitter.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/GlobalFitter.cc $(ROOTFLAG)

#TFitParticleEtEtaPhi.o: TFitParticleEtEtaPhi.cc
#	$(CC) $(CFLAGS) TFitParticleEtEtaPhi.cc $(ROOTFLAG)
#
#TAbsFitParticle.o: TAbsFitParticle.cc
#	$(CC) $(CFLAGS) TAbsFitParticle.cc $(ROOTFLAG)
#
#TFitConstraintM.o: TFitConstraintM.cc
#	$(CC) $(CFLAGS) TFitConstraintM.cc $(ROOTFLAG)
#
#TAbsFitConstraint.o: TAbsFitConstraint.cc 
#	$(CC) $(CFLAGS) TAbsFitConstraint.cc $(ROOTFLAG)
#
#TKinFitter.o: TKinFitter.cc 
#	$(CC) $(CFLAGS) TKinFitter.cc $(ROOTFLAG)

DrawBase.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C fitTools.o
	$(CC) $(CFLAGS) $(INCLUDES) fitTools.o $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C $(ROOTFLAG) $(EXTRALIBS) 

fitTools.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C $(ROOTFLAG) $(EXTRALIBS)

QGLikelihoodCalculator.o: $(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/ $(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C $(ROOTFLAG)

HelicityLikelihoodDiscriminant.o: $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.cc RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/ RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.cc $(ROOTFLAG)

RooBphione.o: $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBphione.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBphione.cc $(ROOTFLAG)

RooBtheta.o: $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBtheta.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBtheta.cc $(ROOTFLAG)

RooPentaSpinTwo.o: $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooPentaSpinTwo.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooPentaSpinTwo.cc $(ROOTFLAG)

RooBkg2l2jHL.o: $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBkg2l2jHL.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBkg2l2jHL.cc $(ROOTFLAG)

RooSpinZero5DV2.o: $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooSpinZero5DV2.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooSpinZero5DV2.cc $(ROOTFLAG)

RooBkgd2L2JV2.o: $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBkgd2L2JV2.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBkgd2L2JV2.cc $(ROOTFLAG)

RooFermi.o: RooFermi.cc
	$(CC) $(CFLAGS) $(INCLUDES) RooFermi.cc $(ROOTFLAG) $(EXTRALIBS)

RooCruijffPdf.o: RooCruijffPdf.cc
	$(CC) $(CFLAGS) $(INCLUDES) RooCruijffPdf.cc $(ROOTFLAG) $(EXTRALIBS)

BTagSFUtil.o: $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/BTagSFUtil.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/ $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/BTagSFUtil.cc $(ROOTFLAG)

PUWeight.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/PUWeight.C
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/PUWeight.C $(ROOTFLAG)


AnalysisJet.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisJet.cc 
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisJet.cc $(ROOTFLAG)

AnalysisElectron.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisElectron.cc AnalysisLepton.o
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisElectron.cc $(ROOTFLAG)

AnalysisMuon.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisMuon.cc AnalysisLepton.o
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisMuon.cc $(ROOTFLAG)

AnalysisLepton.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisLepton.cc
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisLepton.cc $(ROOTFLAG)

AnalysisNeutrino.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisNeutrino.h
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisNeutrino.h $(ROOTFLAG)



clean:
	rm *.o
