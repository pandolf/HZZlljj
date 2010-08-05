CC = g++
CFLAGS = -Wall -c
ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`


do2ndLevel_DiJet: AnalysisJet.o TreeAnalyzer.o TreeAnalyzer_DiJet.o
	$(CC) -Wall -o do2ndLevel_DiJet do2ndLevel_DiJet.cpp AnalysisJet.o TreeAnalyzer.o TreeAnalyzer_DiJet.o $(ROOTFLAG)

do2ndLevel_HZZlljj: Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o 
	$(CC) -Wall -o do2ndLevel_HZZlljj do2ndLevel_HZZlljj.cpp Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o $(ROOTFLAG)

prova: AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o
	$(CC) -Wall -o prova prova_TreeAnalyzer_Photon.cpp AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o $(ROOTFLAG)


Ntp1Analyzer.o: Ntp1Analyzer.C
	$(CC) $(CFLAGS) Ntp1Analyzer.C $(ROOTFLAG)

Ntp1Analyzer_HZZlljj.o: Ntp1Analyzer_HZZlljj.C
	$(CC) $(CFLAGS) Ntp1Analyzer_HZZlljj.C $(ROOTFLAG)

fitTools.o: fitTools.C
	$(CC) $(CFLAGS) fitTools.C $(ROOTFLAG)


clean:
	rm *.o
