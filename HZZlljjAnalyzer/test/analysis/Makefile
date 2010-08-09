CC = g++
CFLAGS = -Wall -c
ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`


merge_and_setWeights: 
	$(CC) -Wall -o merge_and_setWeights merge_and_setWeights.cpp $(ROOTFLAG)

do2ndLevel_HZZlljj: Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o 
	$(CC) -Wall -o do2ndLevel_HZZlljj do2ndLevel_HZZlljj.cpp Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o $(ROOTFLAG)

drawHZZlljj: DrawBase.o fitTools.o drawHZZlljj.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawHZZlljj drawHZZlljj.cpp DrawBase.o fitTools.o $(ROOTFLAG)



Ntp1Analyzer.o: Ntp1Analyzer.C
	$(CC) $(CFLAGS) Ntp1Analyzer.C $(ROOTFLAG)

Ntp1Analyzer_HZZlljj.o: Ntp1Analyzer_HZZlljj.C
	$(CC) $(CFLAGS) Ntp1Analyzer_HZZlljj.C $(ROOTFLAG)

DrawBase.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C $(ROOTFLAG)

fitTools.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C $(ROOTFLAG)



clean:
	rm *.o
	rm merge_and_setWeights
