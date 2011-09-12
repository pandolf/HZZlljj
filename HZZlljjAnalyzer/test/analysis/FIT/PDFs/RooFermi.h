/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_FERMI
#define ROO_FERMI

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

 
class RooFermi : public RooAbsPdf {
public:
  RooFermi();
  RooFermi(const char *name, const char *title,
 	    RooAbsReal& _x,
            RooAbsReal& _cutOff,
	   RooAbsReal& _beta
	   );
  RooFermi(const RooFermi& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooFermi(*this,newname); }
  inline virtual ~RooFermi() { }

protected:

  RooRealProxy x ;
  RooRealProxy cutOff ;
  RooRealProxy beta ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooFermi,1) 
};
 
#endif
