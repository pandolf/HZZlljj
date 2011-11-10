CC = g++
CFLAGS = -Wall -c -g

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/HiggsAnalysis/CombinedLimit/interface

ROOTSYS  ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs` 

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -L$(CMSSW_BASE)/lib/slc5_amd64_gcc434/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit



fitSidebands: fitSidebands.cpp SidebandFitter.o fitTools.o
	$(CC) -Wall $(INCLUDES) -o fitSidebands fitSidebands.cpp SidebandFitter.o fitTools.o $(ROOTFLAG) $(EXTRALIBS) -lHiggsAnalysisCombinedLimit

create_datacards: create_datacards.cpp SidebandFitter.o fitTools.o
	$(CC) -Wall $(INCLUDES) -o create_datacards create_datacards.cpp SidebandFitter.o fitTools.o $(ROOTFLAG) $(EXTRALIBS) -lHiggsAnalysisCombinedLimit

drawMZZ_with_curve: drawMZZ_with_curve.cpp DrawBase.o fitTools.o SidebandFitter.o 
	$(CC) -Wall $(INCLUDES) -o drawMZZ_with_curve drawMZZ_with_curve.cpp DrawBase.o fitTools.o SidebandFitter.o $(ROOTFLAG) $(EXTRALIBS) -lHiggsAnalysisCombinedLimit

draw_golfcourse: draw_golfcourse.cpp DrawBase.o fitTools.o
	$(CC) -Wall $(INCLUDES) -o draw_golfcourse draw_golfcourse.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS) -lHiggsAnalysisCombinedLimit

do_latexYieldTable: do_latexYieldTable.cpp SidebandFitter.o
	$(CC) -Wall $(INCLUDES) -o do_latexYieldTable do_latexYieldTable.cpp SidebandFitter.o $(ROOTFLAG) $(EXTRALIBS) -lHiggsAnalysisCombinedLimit


DrawBase.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C fitTools.o
	$(CC) $(CFLAGS) $(INCLUDES) fitTools.o $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C $(ROOTFLAG) $(EXTRALIBS) 

fitTools.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C $(ROOTFLAG) $(EXTRALIBS)


SidebandFitter.o: SidebandFitter.cc fitTools.o
	$(CC) $(CFLAGS) $(INCLUDES) SidebandFitter.cc fitTools.o $(ROOTFLAG) $(EXTRALIBS) -lHiggsAnalysisCombinedLimit




clean:
	rm *.o
