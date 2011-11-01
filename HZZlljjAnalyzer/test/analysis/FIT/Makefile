CC = g++
CFLAGS = -Wall -c -g

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/

ROOTSYS  ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs` 

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit



fitSidebands: fitSidebands.cpp RooFermi.o RooCB.o SidebandFitter.o PdfDiagonalizer.o
	$(CC) -Wall $(INCLUDES) -o fitSidebands fitSidebands.cpp RooFermi.o RooCB.o SidebandFitter.o PdfDiagonalizer.o $(ROOTFLAG) $(EXTRALIBS)

generateToysForParErrors: generateToysForParErrors.cpp RooFermi.o RooCB.o SidebandFitter.o
	$(CC) -Wall $(INCLUDES) -o generateToysForParErrors generateToysForParErrors.cpp RooFermi.o RooCB.o SidebandFitter.o $(ROOTFLAG) $(EXTRALIBS)

create_datacards: create_datacards.cpp RooFermi.o RooCB.o RooDoubleCB.o RooRelBW.o
	$(CC) -Wall $(INCLUDES) -o create_datacards create_datacards.cpp RooFermi.o RooCB.o RooDoubleCB.o RooRelBW.o $(ROOTFLAG) $(EXTRALIBS)

modifyFitParamErrors: modifyFitParamErrors.cpp SidebandFitter.o RooFermi.o RooCB.o 
	$(CC) -Wall $(INCLUDES) -o modifyFitParamErrors modifyFitParamErrors.cpp SidebandFitter.o RooFermi.o RooCB.o $(ROOTFLAG) $(EXTRALIBS)

drawMZZ_with_curve: drawMZZ_with_curve.cpp DrawBase.o RooFermi.o RooCB.o fitTools.o
	$(CC) -Wall $(INCLUDES) -o drawMZZ_with_curve drawMZZ_with_curve.cpp RooFermi.o RooCB.o DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)


DrawBase.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C fitTools.o
	$(CC) $(CFLAGS) $(INCLUDES) fitTools.o $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C $(ROOTFLAG) $(EXTRALIBS) 

fitTools.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C $(ROOTFLAG) $(EXTRALIBS)


SidebandFitter.o: SidebandFitter.cc
	$(CC) $(CFLAGS) $(INCLUDES) SidebandFitter.cc $(ROOTFLAG) $(EXTRALIBS)

PdfDiagonalizer.o: $(CMSSW_BASE)/src/HiggsAnalysis/CombinedLimit/src/PdfDiagonalizer.cc 
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/HiggsAnalysis/CombinedLimit/src/PdfDiagonalizer.cc $(ROOTFLAG) $(EXTRALIBS)


RooFermi.o: RooFermi.cc
	$(CC) $(CFLAGS) $(INCLUDES) RooFermi.cc $(ROOTFLAG) $(EXTRALIBS)

RooCB.o: RooCB.cc
	$(CC) $(CFLAGS) $(INCLUDES) RooCB.cc $(ROOTFLAG) $(EXTRALIBS)

RooDoubleCB.o: RooDoubleCB.cc
	$(CC) $(CFLAGS) $(INCLUDES) RooDoubleCB.cc $(ROOTFLAG) $(EXTRALIBS)

RooRelBW.o: RooRelBW.cc
	$(CC) $(CFLAGS) $(INCLUDES) RooRelBW.cc $(ROOTFLAG) $(EXTRALIBS)




clean:
	rm *.o
