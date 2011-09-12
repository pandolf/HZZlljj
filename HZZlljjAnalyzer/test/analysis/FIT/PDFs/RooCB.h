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

#ifndef ROO_CB
#define ROO_CB

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;

 
class RooCB : public RooAbsPdf {
public:
  RooCB();
  RooCB(const char *name, const char *title,
	RooAbsReal& _x,
	RooAbsReal& _mean,
	RooAbsReal& _width,
	RooAbsReal& _alpha,
	RooAbsReal& _n,
	RooAbsReal& _theta
	   );
  RooCB(const RooCB& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooCB(*this,newname); }
  inline virtual ~RooCB() { }

protected:

  RooRealProxy x ;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha;
  RooRealProxy n;
  RooRealProxy theta;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooCB,1)
};
 
#endif
