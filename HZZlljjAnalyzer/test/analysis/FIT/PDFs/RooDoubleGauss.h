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

#ifndef ROO_DOUBLEGAUSS
#define ROO_DOUBLEGAUSS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;
 
class RooDoubleGauss : public RooAbsPdf {
public:
  RooDoubleGauss(const char *name, const char *title,
		 RooAbsReal& _x,
		 RooAbsReal& _mean,
		 RooAbsReal& _shift,
		 RooAbsReal& _width,
		 RooAbsReal& _alpha,
		 RooAbsReal& _f

	   );
  RooDoubleGauss(const RooDoubleGauss& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooDoubleGauss(*this,newname); }
  inline virtual ~RooDoubleGauss() { }

protected:

  RooRealProxy x ;
  RooRealProxy mean ;
  RooRealProxy shift ;
  RooRealProxy width ;
  RooRealProxy alpha ;
  RooRealProxy f ;
  Double_t evaluate() const ;

private:

  ClassDef(RooDoubleGauss,0) // Your description goes here...
};
 
#endif
