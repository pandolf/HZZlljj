gC = g++
CFLAGS = -Wall -c
ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`


merge_and_setWeights: merge_and_setWeights.cpp
	$(CC) -Wall -o merge_and_setWeights merge_and_setWeights.cpp $(ROOTFLAG)

do2ndLevel_HZZlljj: Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o fitTools.o do2ndLevel_HZZlljj.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_HZZlljj do2ndLevel_HZZlljj.cpp Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o fitTools.o $(ROOTFLAG)

finalize_HZZlljj: Ntp1Finalizer.o Ntp1Finalizer_HZZlljj.o finalize_HZZlljj.cpp fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitParticleEtEtaPhi.o TKinFitter.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o finalize_HZZlljj finalize_HZZlljj.cpp Ntp1Finalizer.o Ntp1Finalizer_HZZlljj.o fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitParticleEtEtaPhi.o TKinFitter.o $(ROOTFLAG)

finalize_ZGamma: Ntp1Finalizer.o Ntp1Finalizer_ZGamma.o finalize_ZGamma.cpp fitTools.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o finalize_ZGamma finalize_ZGamma.cpp Ntp1Finalizer.o Ntp1Finalizer_ZGamma.o fitTools.o $(ROOTFLAG)

do2ndLevel_JetStudies: Ntp1Analyzer.o Ntp1Analyzer_JetStudies.o 
	$(CC) -Wall -o do2ndLevel_JetStudies do2ndLevel_JetStudies.cpp Ntp1Analyzer.o Ntp1Analyzer_JetStudies.o $(ROOTFLAG)

do2ndLevel_ZGamma: Ntp1Analyzer.o Ntp1Analyzer_ZGamma.o 
	$(CC) -Wall -o do2ndLevel_ZGamma do2ndLevel_ZGamma.cpp Ntp1Analyzer.o Ntp1Analyzer_ZGamma.o $(ROOTFLAG)

do2ndLevel_TMVA: Ntp1Analyzer.o Ntp1Analyzer_TMVA.o do2ndLevel_TMVA.cpp
	$(CC) -Wall -o do2ndLevel_TMVA do2ndLevel_TMVA.cpp Ntp1Analyzer.o Ntp1Analyzer_TMVA.o $(ROOTFLAG)

drawHZZlljj: DrawBase.o fitTools.o drawHZZlljj.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawHZZlljj drawHZZlljj.cpp DrawBase.o fitTools.o $(ROOTFLAG)

drawHZZlljj_kinfit: DrawBase.o fitTools.o drawHZZlljj_kinfit.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawHZZlljj_kinfit drawHZZlljj_kinfit.cpp DrawBase.o fitTools.o $(ROOTFLAG)

drawHZZlljj_signalOnly: DrawBase.o fitTools.o drawHZZlljj_signalOnly.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawHZZlljj_signalOnly drawHZZlljj_signalOnly.cpp DrawBase.o fitTools.o $(ROOTFLAG)

drawHZZlljj_Run2010B: DrawBase.o fitTools.o drawHZZlljj_Run2010B.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawHZZlljj_Run2010B drawHZZlljj_Run2010B.cpp DrawBase.o fitTools.o $(ROOTFLAG)

drawMCHistos: DrawBase.o fitTools.o drawMCHistos.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawMCHistos drawMCHistos.cpp DrawBase.o fitTools.o $(ROOTFLAG)

drawJetStudies: DrawBase.o fitTools.o drawJetStudies.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawJetStudies drawJetStudies.cpp DrawBase.o fitTools.o $(ROOTFLAG)

drawZGamma: DrawBase.o fitTools.o drawZGamma.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawZGamma drawZGamma.cpp DrawBase.o fitTools.o $(ROOTFLAG)

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

Ntp1Finalizer.o: Ntp1Finalizer.C
	$(CC) $(CFLAGS) Ntp1Finalizer.C $(ROOTFLAG)

Ntp1Finalizer_HZZlljj.o: Ntp1Finalizer_HZZlljj.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Finalizer_HZZlljj.C $(ROOTFLAG)

Ntp1Finalizer_ZGamma.o: Ntp1Finalizer_ZGamma.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Finalizer_ZGamma.C $(ROOTFLAG)

Ntp1Finalizer_HZZlljj_TEST.o: Ntp1Finalizer_HZZlljj_TEST.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools Ntp1Finalizer_HZZlljj_TEST.C $(ROOTFLAG)

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



clean:
	rm *.o
