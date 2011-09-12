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

#ifndef ROO_RELBW
#define ROO_RELBW

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;

 
class RooRelBW : public RooAbsPdf {
public:
  RooRelBW();
  RooRelBW(const char *name, const char *title,
	   RooAbsReal& _x,
	   RooAbsReal& _mean,
	   RooAbsReal& _width,
	   RooAbsReal& _n
	   );
  RooRelBW(const RooRelBW& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooRelBW(*this,newname); }
  inline virtual ~RooRelBW() { }

protected:

  RooRealProxy x ;
  RooRealProxy mean ;
  RooRealProxy width ;
  RooRealProxy n ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooRelBW,1)
};
 
#endif
