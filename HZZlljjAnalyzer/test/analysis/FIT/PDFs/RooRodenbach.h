#ifndef ROO_RODENBACH
#define ROO_RODENBACH

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;

 
class RooRodenbach : public RooAbsPdf {
public:
  RooRodenbach();
  RooRodenbach(const char *name, const char *title,
	       RooAbsReal& _x,
	       RooAbsReal& _mean,
	       RooAbsReal& _width,
	       RooAbsReal& _alpha,
	       RooAbsReal& _theta
	   );
  RooRodenbach(const RooRodenbach& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooRodenbach(*this,newname); }
  inline virtual ~RooRodenbach() { }

protected:

  RooRealProxy x ;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha;
  RooRealProxy theta;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooRodenbach,1) 
};
 
#endif
