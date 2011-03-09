CC = g++
CFLAGS = -Wall -c
ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`


merge_and_setWeights: merge_and_setWeights.cpp
	$(CC) -Wall -o merge_and_setWeights merge_and_setWeights.cpp $(ROOTFLAG)

do2ndLevel_HZZlljj: Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o fitTools.o do2ndLevel_HZZlljj.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_HZZlljj do2ndLevel_HZZlljj.cpp Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o fitTools.o $(ROOTFLAG)

finalize_HZZlljj: Ntp1Finalizer.o Ntp1Finalizer_HZZlljj.o finalize_HZZlljj.cpp fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitParticleEtEtaPhi.o TKinFitter.o QGLikelihoodCalculator.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood -o finalize_HZZlljj finalize_HZZlljj.cpp Ntp1Finalizer.o Ntp1Finalizer_HZZlljj.o fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitParticleEtEtaPhi.o TKinFitter.o QGLikelihoodCalculator.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o $(ROOTFLAG) -lMinuit

finalize_ZGamma: Ntp1Finalizer.o Ntp1Finalizer_ZGamma.o finalize_ZGamma.cpp fitTools.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o finalize_ZGamma finalize_ZGamma.cpp Ntp1Finalizer.o Ntp1Finalizer_ZGamma.o fitTools.o $(ROOTFLAG)

finalize_JESsyst: Ntp1Finalizer.o Ntp1Finalizer_JESsyst.o finalize_JESsyst.cpp fitTools.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o finalize_JESsyst finalize_JESsyst.cpp Ntp1Finalizer.o Ntp1Finalizer_JESsyst.o fitTools.o $(ROOTFLAG)

finalize_QG: Ntp1Finalizer.o Ntp1Finalizer_QG.o finalize_QG.cpp fitTools.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o finalize_QG finalize_QG.cpp Ntp1Finalizer.o Ntp1Finalizer_QG.o fitTools.o $(ROOTFLAG) -lMinuit

finalize_QGTMVA: Ntp1Finalizer.o Ntp1Finalizer_QGTMVA.o finalize_QGTMVA.cpp fitTools.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o finalize_QGTMVA finalize_QGTMVA.cpp Ntp1Finalizer.o Ntp1Finalizer_QGTMVA.o fitTools.o $(ROOTFLAG) -lMinuit

finalize_ComputeQGLikelihood: Ntp1Finalizer.o Ntp1Finalizer_ComputeQGLikelihood.o QGLikelihoodCalculator.o finalize_ComputeQGLikelihood.cpp fitTools.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf -o finalize_ComputeQGLikelihood finalize_ComputeQGLikelihood.cpp Ntp1Finalizer.o Ntp1Finalizer_ComputeQGLikelihood.o QGLikelihoodCalculator.o fitTools.o $(ROOTFLAG) -lMinuit

finalize_GG: Ntp1Finalizer.o Ntp1Finalizer_GG.o finalize_GG.cpp fitTools.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o finalize_GG finalize_GG.cpp Ntp1Finalizer.o Ntp1Finalizer_GG.o fitTools.o $(ROOTFLAG) -lMinuit

do2ndLevel_JetStudies: Ntp1Analyzer.o Ntp1Analyzer_JetStudies.o 
	$(CC) -Wall -o do2ndLevel_JetStudies do2ndLevel_JetStudies.cpp Ntp1Analyzer.o Ntp1Analyzer_JetStudies.o $(ROOTFLAG)

do2ndLevel_ZGamma: Ntp1Analyzer.o Ntp1Analyzer_ZGamma.o 
	$(CC) -Wall -o do2ndLevel_ZGamma do2ndLevel_ZGamma.cpp Ntp1Analyzer.o Ntp1Analyzer_ZGamma.o $(ROOTFLAG)

do2ndLevel_TMVA: Ntp1Analyzer.o Ntp1Analyzer_TMVA.o do2ndLevel_TMVA.cpp
	$(CC) -Wall -o do2ndLevel_TMVA do2ndLevel_TMVA.cpp Ntp1Analyzer.o Ntp1Analyzer_TMVA.o $(ROOTFLAG)

do2ndLevel_QG: Ntp1Analyzer.o Ntp1Analyzer_QG.o do2ndLevel_QG.cpp
	$(CC) -Wall -o do2ndLevel_QG do2ndLevel_QG.cpp Ntp1Analyzer.o Ntp1Analyzer_QG.o $(ROOTFLAG)

drawHZZlljj: DrawBase.o fitTools.o drawHZZlljj.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawHZZlljj drawHZZlljj.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawHZZlljj_DATA: DrawBase.o fitTools.o drawHZZlljj_DATA.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawHZZlljj_DATA drawHZZlljj_DATA.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawGG_DATA: DrawBase.o fitTools.o drawGG_DATA.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawGG_DATA drawGG_DATA.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawMCHistos: DrawBase.o fitTools.o drawMCHistos.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawMCHistos drawMCHistos.cpp DrawBase.o fitTools.o $(ROOTFLAG)

drawJetStudies: DrawBase.o fitTools.o drawJetStudies.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawJetStudies drawJetStudies.cpp DrawBase.o fitTools.o $(ROOTFLAG)

drawZGamma: DrawBase.o fitTools.o drawZGamma.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawZGamma drawZGamma.cpp DrawBase.o fitTools.o $(ROOTFLAG)

drawZGamma_DATA: DrawBase.o fitTools.o drawZGamma_DATA.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawZGamma_DATA drawZGamma_DATA.cpp DrawBase.o fitTools.o $(ROOTFLAG)

drawQG: DrawBase.o fitTools.o drawQG.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawQG drawQG.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawRej_vs_Eff: DrawBase.o fitTools.o drawRej_vs_Eff.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawRej_vs_Eff drawRej_vs_Eff.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawRej_vs_Eff_LD: DrawBase.o fitTools.o drawRej_vs_Eff_LD.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawRej_vs_Eff_LD drawRej_vs_Eff_LD.cpp DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawJetStudies_H: DrawBase.o fitTools.o drawJetStudies_H.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawJetStudies_H drawJetStudies_H.cpp DrawBase.o fitTools.o $(ROOTFLAG)



Ntp1Analyzer.o: Ntp1Analyzer.C
	$(CC) $(CFLAGS) Ntp1Analyzer.C $(ROOTFLAG)

Ntp1Analyzer_HZZlljj.o: Ntp1Analyzer_HZZlljj.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Analyzer_HZZlljj.C $(ROOTFLAG)

Ntp1Analyzer_JetStudies.o: Ntp1Analyzer_JetStudies.C
	$(CC) $(CFLAGS) Ntp1Analyzer_JetStudies.C $(ROOTFLAG)

Ntp1Analyzer_ZGamma.o: Ntp1Analyzer_ZGamma.C
	$(CC) $(CFLAGS) Ntp1Analyzer_ZGamma.C $(ROOTFLAG)

Ntp1Analyzer_TMVA.o: Ntp1Analyzer_TMVA.C
	$(CC) $(CFLAGS) Ntp1Analyzer_TMVA.C $(ROOTFLAG)

Ntp1Analyzer_QG.o: Ntp1Analyzer_QG.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Analyzer_QG.C $(ROOTFLAG)

Ntp1Finalizer.o: Ntp1Finalizer.C
	$(CC) $(CFLAGS) Ntp1Finalizer.C $(ROOTFLAG)

Ntp1Finalizer_HZZlljj.o: Ntp1Finalizer_HZZlljj.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood Ntp1Finalizer_HZZlljj.C $(ROOTFLAG)

Ntp1Finalizer_ZGamma.o: Ntp1Finalizer_ZGamma.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Finalizer_ZGamma.C $(ROOTFLAG)

Ntp1Finalizer_QG.o: Ntp1Finalizer_QG.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Finalizer_QG.C $(ROOTFLAG)

Ntp1Finalizer_QGTMVA.o: Ntp1Finalizer_QGTMVA.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Finalizer_QGTMVA.C $(ROOTFLAG)

Ntp1Finalizer_ComputeQGLikelihood.o: Ntp1Finalizer_ComputeQGLikelihood.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf Ntp1Finalizer_ComputeQGLikelihood.C $(ROOTFLAG)

Ntp1Finalizer_GG.o: Ntp1Finalizer_GG.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Finalizer_GG.C $(ROOTFLAG)

Ntp1Finalizer_JESsyst.o: Ntp1Finalizer_JESsyst.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Finalizer_JESsyst.C $(ROOTFLAG)

TFitParticleEtEtaPhi.o: TFitParticleEtEtaPhi.cc
	$(CC) $(CFLAGS) TFitParticleEtEtaPhi.cc $(ROOTFLAG)

TAbsFitParticle.o: TAbsFitParticle.cc
	$(CC) $(CFLAGS) TAbsFitParticle.cc $(ROOTFLAG)

TFitConstraintM.o: TFitConstraintM.cc
	$(CC) $(CFLAGS) TFitConstraintM.cc $(ROOTFLAG)

TAbsFitConstraint.o: TAbsFitConstraint.cc 
	$(CC) $(CFLAGS) TAbsFitConstraint.cc $(ROOTFLAG)

TKinFitter.o: TKinFitter.cc 
	$(CC) $(CFLAGS) TKinFitter.cc $(ROOTFLAG)

DrawBase.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C $(ROOTFLAG)

fitTools.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C $(ROOTFLAG)

QGLikelihoodCalculator.o: $(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/ $(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C $(ROOTFLAG)

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




clean:
	rm *.o
